#%%
from classix import CLASSIX
import numpy as np
from time import time
from tqdm import tqdm
from copy import deepcopy
from collections import deque
import scipy.sparse as sparse
from scipy import sparse
from fractions import Fraction
from spmv import spsubmatxvec

class CLASSIX_T:
    def __init__(self, sorting="popcount", radius=0.3, minPts=1, group_merging="tanimoto_distance", norm=True, mergeScale=1.4, post_alloc=True, mergeTinyGroups=True, verbose=1, short_log_form=True, use_spmv=True):
        self.group_merging = group_merging
        self.mergeScale = mergeScale
        self.mergeTinyGroups = mergeTinyGroups
        self.sorting = sorting
        self.radius = radius
        self.minPts = minPts
        self.norm = norm

        self.post_alloc = post_alloc
        self.mergeTinyGroups = mergeTinyGroups
        self.verbose = verbose
        self.truncate = short_log_form


    
    def __str__(self):
        return f"CLASSIX_T(sorting={self.sorting}, radius={self.radius}, minPts={self.minPts}, group_merging={self.group_merging}, mergeScale={self.mergeScale}, mergeTinyGroups={self.mergeTinyGroups})"
    
    def __repr__(self):
        return f"CLASSIX_T(sorting={self.sorting}, radius={self.radius}, minPts={self.minPts}, group_merging={self.group_merging}, mergeScale={self.mergeScale}, mergeTinyGroups={self.mergeTinyGroups})"
    
    def fit(self, data, r=None, mergeScale=None, minPts=None):
        total_time = time()
        if r is not None:
            self.radius = r
        
        if mergeScale is not None:
            self.mergeScale = mergeScale

        if minPts is not None:
            self.minPts = minPts

        n, fdim = data.shape
        sort_vals = np.sum(data, axis=1)

        # print(sort_vals)

        self.ind = np.argsort(sort_vals, kind='stable')
        self.unsort_ind = np.argsort(self.ind, kind='stable')
        sort_vals = sort_vals[self.ind] 
        data = data[self.ind,:] # sort data

        # print(data)

        print("\nOWN AGGREGATION")


        ips_time = 0.0
        search_time = 0.0
        loop_time = 0.0
        conversion_time = 0.0
        
        st = time()
        lab = 0

        # labels of the data points. Initially all are -1 to keep track of unallocated points
        self.labels = np.array([-1]*n)
        nr_dist = 0
        # Indices of group starting pts, equivalent of the gc list in the original code
        self.splist = []
        rhs = 1/(1-self.radius)+1; # rhs of tanimoto distance inequality
        rhsi = 1/rhs

        rhs_n = 2 - self.radius
        rhs_d = 1 - self.radius

        # Do a small check to see if radius is a rational number:
        # If it is, we can use the rhs_n and rhs_d values to avoid floating point errors!

        self.group_sizes = [] # group sizes
        self.group_starting_pts = [] # group starting pts in the sorted array

        datas = sparse.csr_matrix(data)
        datas_data = datas.data
        datas_indices = datas.indices
        datas_indptr = datas.indptr
        self.aggregation_dict = {}
        
        
        # Aggregation
        for i in tqdm(range(n)):
            
            if self.labels[i] >= 0:
                continue
                
            clustc = data[i,:]
            
            # What if we take the starting point as a slice of the sparse array instead?
            # clustc = datas[i,:]
            #
            self.labels[i] = lab

            # num_group = 1
            self.splist.append(i)
            self.group_sizes.append(1)
            
            search_radius = sort_vals[i]/(1-self.radius) # Only the right side limit is needed
            
            st1 = time()
            last_j = np.searchsorted(sort_vals, search_radius, side='right')
            self.aggregation_dict[i] = last_j

            search_time += time() - st1

            st1 = time()

            # ips = data[i+1:last_j,:]@clustc.T
            ips = spsubmatxvec(datas_data, datas_indptr, datas_indices, i+1, last_j, clustc)

            ips_time += time() - st1
            # Let's make the product with rhsi times sort_vals the condition.
            # lhs = sort_vals[i+1:last_j]
            
            nr_dist += last_j - i - 1
            
            st_c = time()

            vec = ips/(sort_vals[i+1:last_j]+sort_vals[i])
            conversion_time += time() - st_c
            
            st1 = time()

            """
            Create a custom C function that takes in a range of indices, a vector, and a scalar value.
            If the value of the vector at index i is greater than 0, it will check if the value is greater than the scalar value,
            
            Maybe just try boolean masks first?
            Note: self.labels has size = n
                  vec has size = last_j - i - 1

            The boolean mask of vec needs to padded with False values to match the size of self.labels

            Then we can take the intersection of the two boolean masks.

            ```
            vec_mask = vec >= rhsi
            assigned_mask = self.labels >= 0
            
            ```

            """

            notAssigned = self.labels < 0
            vec_mask = vec >= rhsi
            vec_mask = np.pad(vec_mask, (i+1, n - last_j), 'constant', constant_values=False)
            reassignMask = np.logical_and(notAssigned, vec_mask)
            self.group_sizes[-1] += np.sum(reassignMask)
            self.labels[reassignMask] = lab


            # for j in range(i+1, last_j):
            #     if self.labels[j] >= 0:
            #         continue

            #     # if sort_vals[j]+sort_vals[i]<=rhs*ips[j-i-1]:
            #     if (sort_vals[j]+sort_vals[i])*rhs_d<=rhs_n*ips[j-i-1]: # No division!
            #         # num_group += 1
            #         self.labels[j] = lab
            #         self.group_sizes[-1] += 1
                            
            loop_time += time() - st1

            lab += 1
            
        print("time for ips:", ips_time)
        print("time for search:", search_time)
        print("time for conversion:", conversion_time)
        print("time for loop:", loop_time)
        print("nr_dist:", nr_dist)

        self.aggregation_labels = deepcopy(np.array(self.labels))
        self.search_time = search_time
        self.ips_time = ips_time
        self.aggregation_time = time()-st
        
        # merging
        st = time()

        self.spdata = data[self.splist,:]
        spdatas = sparse.csr_matrix(self.spdata)
        spdatas_data = spdatas.data
        spdatas_indices = spdatas.indices
        spdatas_indptr = spdatas.indptr
        spdata_len = len(self.splist)
        self.group_centers = self.spdata

        # Scores of the group starting points
        sort_vals_sp = sort_vals[self.splist]

        # Group labels of the group starting points
        label_sp = self.labels[self.splist]
        scale = self.mergeScale

        
        

        self.Adj = np.zeros((len(self.splist), len(self.splist)), dtype=np.int8)
        
        t1 = 0.0
        t2 = 0.0
        t3 = 0.0
        minPts = self.minPts

        for i in tqdm(range(len(self.splist))):
            if not self.mergeTinyGroups:
                if self.group_sizes[i] < minPts:
                    continue
                
            xi = self.spdata[i, :]
            search_radius = sort_vals_sp[i]/(1-self.mergeScale*self.radius)

            st1 = time()
            last_j = np.searchsorted(sort_vals_sp, search_radius, side='right')
            first_j = i
            ips = spsubmatxvec(spdatas_data, spdatas_indptr, spdatas_indices, first_j, last_j, xi)
            # ips = np.matmul(self.spdata[first_j:last_j,:],xi.T)
            tanimoto_distance = 1 - ips/(sort_vals_sp[i] + sort_vals_sp[first_j:last_j] - ips)
            # tanimoto_distance = 1 - tanimoto_score
            t1 += time() - st1

            st2 = time()
            inds = np.where(tanimoto_distance <= scale*self.radius)

            if not self.mergeTinyGroups:
                inds = inds[0][self.group_sizes[inds[0]] >= minPts]
                inds = i + inds
            
            else:
                inds = i + inds[0]



            self.Adj[i, inds] = 1
            self.Adj[inds, i] = 1

            spl = np.unique(label_sp[ inds ])

            minlab = np.min(spl)

            t2 += time() - st2

            st3 = time()
            for label in spl:
                label_sp[label_sp==label] = minlab
            t3 += time() - st3
            # continue

        self.initial_merging_labels = deepcopy(label_sp)

        print("  merging time:", time()-st)
        self.merging_time = time()-st
        
        # Cluster redistribution
        print(" minPts Merging")
        st = time()
        ############################################################################################################

        labels_sp_copy = deepcopy(label_sp)
        ul = np.unique(label_sp)
        cs = np.zeros(len(ul))
        self.group_sizes = np.array(self.group_sizes)
        for i in range(len(ul)):
            id = np.where(label_sp==ul[i])
            label_sp[id] = i
            cs[i] = np.sum(self.group_sizes[id])

        small_clusters = np.where(cs<self.minPts)[0]

        # print("small clusters", small_clusters)

        labels_sp_copy_2 = deepcopy(label_sp)

        for i in small_clusters:
            ii = np.where(labels_sp_copy_2==i)[0]
            for iii in ii:
                xi = self.spdata[iii, :]
                # d = Tanimoto distances to all other group centers
                # ips = np.matmul(self.spdata, xi.T)
                ips = spsubmatxvec(spdatas_data, spdatas_indptr, spdatas_indices, 0, spdata_len, xi)
                d = 1 - ips/(sort_vals_sp[iii] + sort_vals_sp - ips)
                
                o = np.argsort(d, kind='stable')
                # merge with the closest group that has more than minPts
                for j in o:
                    if cs[labels_sp_copy_2[j]]>=minPts:
                        label_sp[iii] = labels_sp_copy_2[j]
                        # Instead of having a separate Adjanceny matrix for the points merged with minPts criteria
                        # we can just update the old Adjacency matrix here. The issue is that the distance matrix will
                        # have some weird discrepancies.
                        self.Adj[iii, j] = 2
                        self.Adj[j, iii] = 2
                        break

        ul = np.unique(label_sp)
        cs = np.zeros(len(ul))
        for i in range(len(ul)):
            id = np.where(label_sp==ul[i])
            label_sp[id] = i
            cs[i] = np.sum(self.group_sizes[id])

        # print("final cluster sizes", cs)

        print("  minPts merging time:", time()-st)
        
        for idx, label in enumerate(self.labels):
            self.labels[idx] = label_sp[label]

        ############################################################################################################

        self.labels = self.labels[self.unsort_ind]
        # print(self.labels)
        self.aggregation_labels = self.aggregation_labels[self.unsort_ind]
        self.group_labels = label_sp
        self.group_centre_pts = self.spdata
        self.group_centers = self.splist
        print("Total time:", time()-total_time)
        return self

    def explain(self, ind1=None, ind2=None):
        if ind1 is None and ind2 is None:
            # If there are no specific points to explain, print the general information
            print(f"The data was clustered into {len(self.group_labels)} groups. These were further merged into {len(np.unique(self.group_labels))} clusters.")

        if ind1 is not None and ind2 is None:
            # If only one index is provided, print the cluster information for that index
            print(f"The data point at index {ind1} was assigned to cluster {self.labels[ind1]}")

        if ind1 is not None and ind2 is not None:
            # If two indices are provided, print the cluster information
            if self.labels[ind1] == self.labels[ind2]:
                # If the two data points belong to the same cluster, check if they are in the same group
                print(f"The data points at indices {ind1} and {ind2} belong to the same cluster.")
                agg_label_1 = self.aggregation_labels[ind1]
                agg_label_2 = self.aggregation_labels[ind2]
                if agg_label_1 == agg_label_2:
                    # If the two data points are in the same group, print the information
                    print(f"The data points at indices {ind1} and {ind2} are in the same group {agg_label_1}.")
                    return None

                else:
                    # If the two data points are in different groups, find the shortest connection path between the two groups
                    print(f"The data points at indices {ind1} and {ind2} are in different groups.")
                    connected_path = self.bfs_shortest_path(self.Adj, agg_label_1, agg_label_2)
                    if connected_path is not None:
                        # If a connection path is found, print the path
                        print(f'The connections between {ind1} and {ind2} are via this path: {connected_path[0]} ', end="")
                        for k in range(len(connected_path)-1):
                            if self.Adj[connected_path[k]-1, connected_path[k+1]-1] == 2:
                                print(f'(minPts) -> {connected_path[k+1]}', end="")
                            else:
                                print(f' -> {connected_path[k+1]}', end="")

                        return connected_path
                    
                    else:
                        # If no connection path is found, print that there is no connection, there must be something wrong with the code
                        print(f"No connection path found between {ind1} and {ind2}, although they belong to different groups in the same cluster. Please check the program for bugs!")
                        return None

            else:
                print(f"The data points at indices {ind1} and {ind2} belong to different clusters.")


    def bfs_shortest_path(self, adj_matrix, start, goal):
        # Convert start and goal from 1-based to 0-based indices
        start -= 1
        goal -= 1
        
        # Initialize the queue with the start node and a path containing only the start node
        queue = deque([(start, [start])])
        
        # Keep track of visited nodes to avoid cycles
        visited = set()
        
        while queue:
            # Dequeue a node and the path that led to it
            current_node, path = queue.popleft()
            
            # If the current node is the goal, return the path
            if current_node == goal:
                return [node + 1 for node in path]  # Convert back to 1-based indices
            
            # Mark the current node as visited
            visited.add(current_node)
            
            # Enqueue all unvisited neighbors
            for neighbor, connected in enumerate(adj_matrix[current_node]):
                if connected and neighbor not in visited:
                    queue.append((neighbor, path + [neighbor]))

        return None
# %%
