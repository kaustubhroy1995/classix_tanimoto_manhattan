void multiply(
    const int *m_data, const int *m_indptr, const int *m_indices,
    int start_row, int end_row, const int *v_data, int *result) {
    
    int j = 0;
    int val = 0;
    int num_vals_matrix = 0;
    int idx = 0;
    int m_indptri = 0;
    
    for (int i = start_row; i < end_row; i++) {
        m_indptri = m_indptr[i];
        num_vals_matrix = m_indptr[i+1] - m_indptri;
        val = 0;
        for (j = m_indptri; j < m_indptri + num_vals_matrix; j++) { 
            idx = m_indices[j];
            val += m_data[j] * v_data[idx];
        }
        result[i-start_row] = val;
    }
}