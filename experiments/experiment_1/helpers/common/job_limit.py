"""Windows job-object commit limit for a paused benchmark process tree.

The controller must start its child with CREATE_NO_WINDOW and keep the actual
worker paused on a READY/GO handshake until assign_job returns. Pass the actual
worker PID from READY, not the outer Windows venv-launcher PID: an already-running
launcher can use a separate nested job for its child, preventing both processes
from being added to a common new job. Only that explicitly identified worker is
assigned; subsequently created children inherit the job. The worker must have
no pre-existing subprocesses that require bounding. Process-tree discovery is
deliberately disabled: Windows parent PID snapshots can identify unrelated old
processes after a parent PID is reused. No breakaway or host-wide limit is
enabled. The outer launcher holds no benchmark arrays.

The memory limit concerns committed memory, not working-set size. Assignment
does not retroactively reject allocations made before joining the job; attach
before substantial allocations when a strict whole-run bound is needed.

API:
    job = assign_job(pid=ready["worker_pid"], max_bytes=4 * 1024**3)
    ... send GO, wait for the child ...
    close_job(job)  # also kills any remaining assigned descendants

Alternatively pass an existing process handle as the first positional argument.
An integer positional argument is a HANDLE, never interpreted as a PID.

References:
https://learn.microsoft.com/en-us/windows/win32/api/winnt/ns-winnt-jobobject_extended_limit_information
https://learn.microsoft.com/en-us/windows/win32/api/jobapi2/nf-jobapi2-assignprocesstojobobject
"""
import ctypes
from ctypes import wintypes
import os

if os.name != "nt":
    raise OSError("This benchmark resource limiter requires Windows")

_kernel = ctypes.WinDLL("kernel32", use_last_error=True)
_SIZE_T = ctypes.c_size_t
_ULONG_PTR = ctypes.c_size_t
_PROCESS_RIGHTS = 0x0001 | 0x0100 | 0x1000  # TERMINATE, SET_QUOTA, QUERY_LIMITED
_JOB_INFO_EXTENDED = 9
_JOB_OBJECT_LIMIT_JOB_MEMORY = 0x00000200
_JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE = 0x00002000


class _BasicLimit(ctypes.Structure):
    _fields_ = [
        ("PerProcessUserTimeLimit", ctypes.c_int64),
        ("PerJobUserTimeLimit", ctypes.c_int64),
        ("LimitFlags", wintypes.DWORD),
        ("MinimumWorkingSetSize", _SIZE_T),
        ("MaximumWorkingSetSize", _SIZE_T),
        ("ActiveProcessLimit", wintypes.DWORD),
        ("Affinity", _ULONG_PTR),
        ("PriorityClass", wintypes.DWORD),
        ("SchedulingClass", wintypes.DWORD),
    ]


class _IoCounters(ctypes.Structure):
    _fields_ = [(name, ctypes.c_uint64) for name in (
        "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
        "ReadTransferCount", "WriteTransferCount", "OtherTransferCount",
    )]


class _ExtendedLimit(ctypes.Structure):
    _fields_ = [
        ("BasicLimitInformation", _BasicLimit),
        ("IoInfo", _IoCounters),
        ("ProcessMemoryLimit", _SIZE_T),
        ("JobMemoryLimit", _SIZE_T),
        ("PeakProcessMemoryUsed", _SIZE_T),
        ("PeakJobMemoryUsed", _SIZE_T),
    ]


_kernel.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
_kernel.CreateJobObjectW.restype = wintypes.HANDLE
_kernel.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
_kernel.SetInformationJobObject.restype = wintypes.BOOL
_kernel.QueryInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD, ctypes.c_void_p]
_kernel.QueryInformationJobObject.restype = wintypes.BOOL
_kernel.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
_kernel.AssignProcessToJobObject.restype = wintypes.BOOL
_kernel.IsProcessInJob.argtypes = [wintypes.HANDLE, wintypes.HANDLE, ctypes.POINTER(wintypes.BOOL)]
_kernel.IsProcessInJob.restype = wintypes.BOOL
_kernel.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
_kernel.OpenProcess.restype = wintypes.HANDLE
_kernel.GetProcessId.argtypes = [wintypes.HANDLE]
_kernel.GetProcessId.restype = wintypes.DWORD
_kernel.CloseHandle.argtypes = [wintypes.HANDLE]
_kernel.CloseHandle.restype = wintypes.BOOL


def _failure(operation):
    code = ctypes.get_last_error()
    detail = ctypes.FormatError(code).strip()
    return OSError(code, f"{operation}: {detail}; benchmark must not proceed without its job limit")


def _open_process(pid):
    handle = _kernel.OpenProcess(_PROCESS_RIGHTS, False, pid)
    if not handle:
        raise _failure(f"OpenProcess({pid})")
    return handle


def process_in_job(job, *, pid):
    """Check a READY worker PID is contained before sending GO."""
    handle = _open_process(pid)
    try:
        present = wintypes.BOOL()
        if not _kernel.IsProcessInJob(handle, job, ctypes.byref(present)):
            raise _failure(f"IsProcessInJob({pid})")
        return bool(present.value)
    finally:
        _kernel.CloseHandle(handle)


def assign_job(process_handle=None, max_bytes=4 * 1024**3, *, pid=None,
               include_existing_descendants=False):
    """Bound the explicit worker; its future children inherit the same job.

    Exactly one of process_handle and pid is required. The worker must have no
    existing subprocesses requiring this limit. include_existing_descendants=True
    is rejected before opening/assigning any process; PID-only tree discovery is
    unsafe when Windows reuses PIDs. A failure, including incompatible
    enclosing/nested jobs, raises; partially assigned processes are terminated.
    The caller owns only the returned job handle, not a passed process handle.
    """
    if (process_handle is None) == (pid is None):
        raise ValueError("Pass exactly one process_handle or keyword pid")
    if include_existing_descendants is not False:
        raise ValueError("Existing-descendant discovery is disabled; pass the actual READY worker with no pre-existing subprocesses")
    if not isinstance(max_bytes, int) or not 0 < max_bytes <= ctypes.c_size_t(-1).value:
        raise ValueError("max_bytes must be a positive SIZE_T integer")
    own_process_handle = process_handle is None
    if own_process_handle:
        process_handle = _open_process(pid)
    job = None
    try:
        root_pid = int(_kernel.GetProcessId(process_handle))
        if not root_pid:
            raise _failure("GetProcessId")
        if root_pid == os.getpid():
            raise ValueError("Refusing to assign the controller itself to a kill-on-close job")
        job = _kernel.CreateJobObjectW(None, None)
        if not job:
            raise _failure("CreateJobObjectW")
        limits = _ExtendedLimit()
        limits.BasicLimitInformation.LimitFlags = (
            _JOB_OBJECT_LIMIT_JOB_MEMORY | _JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        )
        limits.JobMemoryLimit = max_bytes
        if not _kernel.SetInformationJobObject(job, _JOB_INFO_EXTENDED, ctypes.byref(limits), ctypes.sizeof(limits)):
            raise _failure("SetInformationJobObject")
        if not _kernel.AssignProcessToJobObject(job, process_handle):
            raise _failure(f"AssignProcessToJobObject({root_pid}); nested-job incompatibility is possible")
        # No enumeration: only the READY process is assigned. Windows ensures
        # subsequently created children inherit the same aggregate job limit.
        return job
    except BaseException:
        if job:
            _kernel.CloseHandle(job)
        raise
    finally:
        if own_process_handle:
            _kernel.CloseHandle(process_handle)


def job_snapshot(job):
    """Return diagnostic Windows job counters.

    Peak counters can include a subsequently denied allocation attempt, so they
    may exceed the configured quota; do not interpret them as resident memory.
    """
    limits = _ExtendedLimit()
    if not _kernel.QueryInformationJobObject(job, _JOB_INFO_EXTENDED, ctypes.byref(limits), ctypes.sizeof(limits), None):
        raise _failure("QueryInformationJobObject")
    return {
        "memory_limit_bytes": int(limits.JobMemoryLimit),
        "peak_job_commit_bytes": int(limits.PeakJobMemoryUsed),
        "peak_process_commit_bytes": int(limits.PeakProcessMemoryUsed),
        "limit_flags": int(limits.BasicLimitInformation.LimitFlags),
    }


def close_job(job):
    """Close a valid job handle; any still-running member processes are killed."""
    if job and not _kernel.CloseHandle(job):
        raise _failure("CloseHandle(job)")
