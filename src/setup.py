from setuptools import setup, dist


class BinaryDistribution(dist.Distribution):
    def is_pure(self):
        return False


setup(
    name='spmv',
    version='0.1',
    description='A faster sparse submatrix - dense vector multiplication in Python',
    package_data={'spmv': ['lib/*.so']},
    include_package_data=True,
    distclass=BinaryDistribution,
    packages=['spmv'],
)