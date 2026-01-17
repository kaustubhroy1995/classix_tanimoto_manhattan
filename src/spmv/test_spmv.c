#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void dot_vector_transposed(
    const int *m_data, const int *m_indptr, const int *m_indices,
    int start_row, int end_row,
    const int *v_data, 
    int *result) {

    for (int i = start_row; i < end_row; i++) {
        int idx_1 = 0;
        int num_vals_matrix = m_indptr[i + 1] - m_indptr[i];
        int val = 0;

        while (idx_1 < num_vals_matrix) {
            int matrix_index = m_indices[m_indptr[i] + idx_1];
            val += m_data[m_indptr[i] + idx_1] * v_data[matrix_index];
            idx_1++;
        }
        result[i] = val;
    }
}
