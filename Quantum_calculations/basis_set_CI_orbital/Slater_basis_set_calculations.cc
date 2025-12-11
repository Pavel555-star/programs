#include "Slater_basis_set_calculations.h"
using namespace std;
template <typename T>
T Slater_basis_set_calculations<T>::Determinant(unsigned int order, T* pointer, T* buffer, T* denominator, T* temp1, T* temp2) 
    { // calculation of determinant matrix (reliable for matrices to 8 * 8 bytes and datetype long long)
    unsigned int i; // pointer, buffer and denominator - order * order dimensions and temp1 and temp2 - order dimensions
    unsigned int j; // for order > 64, else nullptr
    T det = 0;
    T line = 1;
    T matrix[3][5];

    switch (order)
        {
        case 0:
            return(-1); // For incorrectly entered matrices return 0
        case 1:
            return(pointer[0]); // For first order return value of element
        case 2:
            return((pointer[0]*pointer[3])-(pointer[1]*pointer[2])); // For second order return determinant of second-order matrix
        case 3:{ // Calculation by Sarus law
            matrix[0][0] = pointer[0]; // Copying of matrix into auxiliary matrix
            matrix[0][1] = pointer[1];
            matrix[0][2] = pointer[2];
            matrix[1][0] = pointer[3];
            matrix[1][1] = pointer[4];
            matrix[1][2] = pointer[5];
            matrix[2][0] = pointer[6];
            matrix[2][1] = pointer[7];
            matrix[2][2] = pointer[8];
            for (i = 0; i < 3; i++) // Filling of auxiliary lines
                { 
                matrix[i][3] = matrix[i][0];
                matrix[i][4] = matrix[i][1];
                }
            for (j = 0; j < 3; j++) // Addition of plus lines
                {
                for (i =0; i < 3; i++)
                    line = line * matrix[i][i + j];
                    
                det = det + line;
                line = 1;
                }      
            for (j = 0; j < 3; j++) // Subtraction of minus lines
                {
                for (i =0; i < 3; i++)
                    line = line * matrix[2 - i][i + j];
                    
                det = det - line;
                line = 1;
                }
        return det;
        }
    }
    if ((order >= 4) and (order <= 7)) // Calculation by development according to the column
        {
        unsigned int x;
        unsigned int subi;
        unsigned int subj;
        int prefactor;
        T* submatrix;
        
        try {
            submatrix = new T[(order-1) * (order-1)];}
        catch (int)
            {
            determinant_exception_handle = 1;
            }
        for (x = 0; x < order; x++)
            {
            subi = 0;
            for (i = 1; i < order; i++) 
                {
                subj = 0;
                for (j = 0; j < order; j++)
                    {
                    if (j == x)
                        continue;
                    submatrix[subi+((order-1) * subj)] = pointer[i + (order * j)];
                    subj++;
                    }
                subi++;
                }
            if (pointer[order * x] != 0)
                {
                if (x % 2 == 0)
                    prefactor = 1;
                else
                    prefactor = -1;
                det = det + (prefactor * pointer[order * x] * this->Determinant(order-1, submatrix,
                nullptr, nullptr, nullptr,nullptr));
                }
            }
        delete[] submatrix;
        }
    // vectorisation code    
    if (order >= 8 and order <= 64) // Gaussian elimination calculation in triangular form
        {
        unsigned int line_G;
        unsigned int k;
        int sign = 1;
        unsigned int red_order;
        T den;
        T n;
        T d;
        T buffer_s[4096];
        T denominator_s[4096];
        T temp1_s[64];
        T temp2_s[64];
        
        for (i = 0; i < (order * order); i++)
            {
            buffer_s[i] = pointer[i]; // Copying of matrix to modifications
            denominator_s[i] = 1; // Filling of denominator_s matrix by ones
            }
        for (i = 0; i < order; i++)
            {
            for (j = i; j < order; j++)
                {
                if (buffer_s[i + (order * j)] != 0) // Finding of first line with non-zero value in the given column
                {
                line_G = j;
                break;
                }
                if ((j == (order-1)) and (buffer_s[i + (order * j)] == 0)) // Zero control
                  return(0);
                }
            for (j = i; j < order; j++) // Addition and subtraction of lines
                {
                if ((line_G != j) and (buffer_s[i + (order * j)] != 0)) 
                    {
                    n = buffer_s[i + (order * j)];
                    d = buffer_s[i + (order * line_G)];
                    if (d == 0)
                        d = 1;
                        
                    if ((n - (n/d) * d) == 0) // If nominator and denominator_s are divisible without remainder
                        {
                        n = n/d;
                        d = 1;
                            for (k = 0; k < order - 7; k+=8)
                                {
                                buffer_s[k + (order * j)] = buffer_s[k + (order * j)] - (buffer_s[k + (order * line_G)] * n);
                                buffer_s[k + (order * j) + 1] = buffer_s[k + (order * j) + 1] - (buffer_s[k + (order * line_G) + 1] * n);
                                buffer_s[k + (order * j) + 2] = buffer_s[k + (order * j) + 2] - (buffer_s[k + (order * line_G) + 2] * n);
                                buffer_s[k + (order * j) + 3] = buffer_s[k + (order * j) + 3] - (buffer_s[k + (order * line_G) + 3] * n);
                                buffer_s[k + (order * j) + 4] = buffer_s[k + (order * j) + 4] - (buffer_s[k + (order * line_G) + 4] * n);
                                buffer_s[k + (order * j) + 5] = buffer_s[k + (order * j) + 5] - (buffer_s[k + (order * line_G) + 5] * n);
                                buffer_s[k + (order * j) + 6] = buffer_s[k + (order * j) + 6] - (buffer_s[k + (order * line_G) + 6] * n);
                                buffer_s[k + (order * j) + 7] = buffer_s[k + (order * j) + 7] - (buffer_s[k + (order * line_G) + 7] * n);
                                }
                            red_order = order % 8;
                            for (k = order - red_order; k < order; k++)
                            {
                                buffer_s[k + (order * j)] = buffer_s[k + (order * j)] - (buffer_s[k + (order * line_G)] * n);
                            }
                        }
                    else // If nominator and denominator_s are not divisible without remainder
                        {
                        for (k = 0; k < order - 7; k+=8)
                            {
                            buffer_s[k + (order * j)] = buffer_s[k + (order * j)] * d;
                            buffer_s[k + (order * j) + 1] = buffer_s[k + (order * j) + 1] * d;
                            buffer_s[k + (order * j) + 2] = buffer_s[k + (order * j) + 2] * d;
                            buffer_s[k + (order * j) + 3] = buffer_s[k + (order * j) + 3] * d;
                            buffer_s[k + (order * j) + 4] = buffer_s[k + (order * j) + 4] * d;
                            buffer_s[k + (order * j) + 5] = buffer_s[k + (order * j) + 5] * d;
                            buffer_s[k + (order * j) + 6] = buffer_s[k + (order * j) + 6] * d;
                            buffer_s[k + (order * j) + 7] = buffer_s[k + (order * j) + 7] * d;
                            buffer_s[k + (order * j)] = buffer_s[k + (order * j)] - (buffer_s[k + (order * line_G)] * n);
                            buffer_s[k + (order * j) + 1] = buffer_s[k + (order * j) + 1] - (buffer_s[k + (order * line_G) + 1] * n);
                            buffer_s[k + (order * j) + 2] = buffer_s[k + (order * j) + 2] - (buffer_s[k + (order * line_G) + 2] * n);
                            buffer_s[k + (order * j) + 3] = buffer_s[k + (order * j) + 3] - (buffer_s[k + (order * line_G) + 3] * n);
                            buffer_s[k + (order * j) + 4] = buffer_s[k + (order * j) + 4] - (buffer_s[k + (order * line_G) + 4] * n);
                            buffer_s[k + (order * j) + 5] = buffer_s[k + (order * j) + 5] - (buffer_s[k + (order * line_G) + 5] * n);
                            buffer_s[k + (order * j) + 6] = buffer_s[k + (order * j) + 6] - (buffer_s[k + (order * line_G) + 6] * n);
                            buffer_s[k + (order * j) + 7] = buffer_s[k + (order * j) + 7] - (buffer_s[k + (order * line_G) + 7] * n);
                            denominator_s[k + (order * j)] = denominator_s[k + (order * j)] * d;
                            denominator_s[k + (order * j) + 1] = denominator_s[k + (order * j) + 1] * d;
                            denominator_s[k + (order * j) + 2] = denominator_s[k + (order * j) + 2] * d;
                            denominator_s[k + (order * j) + 3] = denominator_s[k + (order * j) + 3] * d;
                            denominator_s[k + (order * j) + 4] = denominator_s[k + (order * j) + 4] * d;
                            denominator_s[k + (order * j) + 5] = denominator_s[k + (order * j) + 5] * d;
                            denominator_s[k + (order * j) + 6] = denominator_s[k + (order * j) + 6] * d;
                            denominator_s[k + (order * j) + 7] = denominator_s[k + (order * j) + 7] * d;
                            }
                        red_order = order % 8;
                        for (k = order - red_order; k < order; k++)
                            {
                            buffer_s[k + (order * j)] = buffer_s[k + (order * j)] * d;
                            buffer_s[k + (order * j)] = buffer_s[k + (order * j)] - (buffer_s[k + (order * line_G)] * n);
                            denominator_s[k + (order * j)] = denominator_s[k + (order * j)] * d;
                            }
                        }
                    }
                }
            if (line_G > i) // Swap of lines
                {
                sign = sign * -1;
                for (k = 0; k < order; k++)
                    {
                    temp1_s[k] = buffer_s[k + (order * line_G)];
                    temp2_s[k] = denominator_s[k + (order * line_G)];
                    buffer_s[k + (order * line_G)] = buffer_s[k + (order * i)];
                    denominator_s[k + (order * line_G)] = denominator_s[k + (order * i)];
                    buffer_s[k + (order * i)] = temp1_s[k];
                    denominator_s[k + (order * i)] = temp2_s[k];
                    }
                }
            }
        det = 1;
        den = 1;
        for (i = 0; i < order; i++)
            {
            det = det * buffer_s[i * (order + 1)];
            den = den * denominator_s[i * (order + 1)];
            if ((det - (det/den) * den) == 0) 
                {
                det = det/den;
                den =1;
                }
            }
        det = (det/den) * sign;
        }
    if (order > 64) // Gaussian elimination calculation in triangular form
        {
        unsigned int line_G;
        unsigned int k;
        int sign = 1;
        unsigned int red_order;
        T den;
        T n;
        T d;
        for (i = 0; i < (order * order); i++)
            {
            buffer[i] = pointer[i]; // Copying of matrix to modifications
            denominator[i] = 1; // Filling of denominator matrix by ones
            }
        for (i = 0; i < order; i++)
            {
            for (j = i; j < order; j++)
                {
                if (buffer[i + (order * j)] != 0) // Finding of first line with non-zero value in the given column
                    {
                    line_G = j;
                    break;
                    }
                if ((j == (order-1)) and (buffer[i + (order * j)] == 0)) // Zero control
                  return(0);
                }
            for (j = i; j < order; j++) // Addition and subtraction of lines
                {
                if ((line_G != j) and (buffer[i + (order * j)] != 0)) 
                    {
                    n = buffer[i + (order * j)];
                    d = buffer[i + (order * line_G)];
                    if (d == 0)
                        d = 1;
                        
                    if ((n - (n/d) * d) == 0) // If nominator and denominator are divisible without remainder
                        {
                        n = n/d;
                        d = 1;
                        for (k = 0; k < order - 7; k+=8)
                            {
                            buffer[k + (order * j)] = buffer[k + (order * j)] - (buffer[k + (order * line_G)] * n);
                            buffer[k + (order * j) + 1] = buffer[k + (order * j) + 1] - (buffer[k + (order * line_G) + 1] * n);
                            buffer[k + (order * j) + 2] = buffer[k + (order * j) + 2] - (buffer[k + (order * line_G) + 2] * n);
                            buffer[k + (order * j) + 3] = buffer[k + (order * j) + 3] - (buffer[k + (order * line_G) + 3] * n);
                            buffer[k + (order * j) + 4] = buffer[k + (order * j) + 4] - (buffer[k + (order * line_G) + 4] * n);
                            buffer[k + (order * j) + 5] = buffer[k + (order * j) + 5] - (buffer[k + (order * line_G) + 5] * n);
                            buffer[k + (order * j) + 6] = buffer[k + (order * j) + 6] - (buffer[k + (order * line_G) + 6] * n);
                            buffer[k + (order * j) + 7] = buffer[k + (order * j) + 7] - (buffer[k + (order * line_G) + 7] * n);
                            }
                        red_order = order % 8;
                        for (k = order - red_order; k < order; k++)
                            buffer[k + (order * j)] = buffer[k + (order * j)] - (buffer[k + (order * line_G)] * n);
                        }
                    else // If nominator and denominator are not divisible without remainder
                        {
                        for (k = 0; k < order - 7; k+=8)
                            {
                            buffer[k + (order * j)] = buffer[k + (order * j)] * d;
                            buffer[k + (order * j) + 1] = buffer[k + (order * j) + 1] * d;
                            buffer[k + (order * j) + 2] = buffer[k + (order * j) + 2] * d;
                            buffer[k + (order * j) + 3] = buffer[k + (order * j) + 3] * d;
                            buffer[k + (order * j) + 4] = buffer[k + (order * j) + 4] * d;
                            buffer[k + (order * j) + 5] = buffer[k + (order * j) + 5] * d;
                            buffer[k + (order * j) + 6] = buffer[k + (order * j) + 6] * d;
                            buffer[k + (order * j) + 7] = buffer[k + (order * j) + 7] * d;
                            buffer[k + (order * j)] = buffer[k + (order * j)] - (buffer[k + (order * line_G)] * n);
                            buffer[k + (order * j) + 1] = buffer[k + (order * j) + 1] - (buffer[k + (order * line_G) + 1] * n);
                            buffer[k + (order * j) + 2] = buffer[k + (order * j) + 2] - (buffer[k + (order * line_G) + 2] * n);
                            buffer[k + (order * j) + 3] = buffer[k + (order * j) + 3] - (buffer[k + (order * line_G) + 3] * n);
                            buffer[k + (order * j) + 4] = buffer[k + (order * j) + 4] - (buffer[k + (order * line_G) + 4] * n);
                            buffer[k + (order * j) + 5] = buffer[k + (order * j) + 5] - (buffer[k + (order * line_G) + 5] * n);
                            buffer[k + (order * j) + 6] = buffer[k + (order * j) + 6] - (buffer[k + (order * line_G) + 6] * n);
                            buffer[k + (order * j) + 7] = buffer[k + (order * j) + 7] - (buffer[k + (order * line_G) + 7] * n);
                            denominator[k + (order * j)] = denominator[k + (order * j)] * d;
                            denominator[k + (order * j) + 1] = denominator[k + (order * j) + 1] * d;
                            denominator[k + (order * j) + 2] = denominator[k + (order * j) + 2] * d;
                            denominator[k + (order * j) + 3] = denominator[k + (order * j) + 3] * d;
                            denominator[k + (order * j) + 4] = denominator[k + (order * j) + 4] * d;
                            denominator[k + (order * j) + 5] = denominator[k + (order * j) + 5] * d;
                            denominator[k + (order * j) + 6] = denominator[k + (order * j) + 6] * d;
                            denominator[k + (order * j) + 7] = denominator[k + (order * j) + 7] * d;
                            }
                        red_order = order % 8;
                        for (k = order - red_order; k < order; k++)
                            {
                            buffer[k + (order * j)] = buffer[k + (order * j)] * d;
                            buffer[k + (order * j)] = buffer[k + (order * j)] - (buffer[k + (order * line_G)] * n);
                            denominator[k + (order * j)] = denominator[k + (order * j)] * d;
                            }
                        }
                    }
                }
            if (line_G > i) // Swap of lines
                {
                sign = sign * -1;
                for (k = 0; k < order; k++)
                    {
                    temp1[k] = buffer[k + (order * line_G)];
                    temp2[k] = denominator[k + (order * line_G)];
                    buffer[k + (order * line_G)] = buffer[k + (order * i)];
                    denominator[k + (order * line_G)] = denominator[k + (order * i)];
                    buffer[k + (order * i)] = temp1[k];
                    denominator[k + (order * i)] = temp2[k];
                    }
                }
            }
        det = 1;
        den = 1;
        for (i = 0; i < order; i++)
            {
            det = det * buffer[i * (order + 1)];
            den = den * denominator[i * (order + 1)];
            if ((det - (det/den) * den) == 0)
                {
                det = det/den;
                den =1;
                }
            }
        det = (det/den)*sign;
        }
    // end of vectorisation code    
    return det;
    }
template <typename T>
int Slater_basis_set_calculations<T>::basis_set_Determinant_set(unsigned int order, T* pointer,unsigned int count, T min, T step,
T* output_values)
    {
    unsigned int i;
    unsigned int j;
    T* diagonal = nullptr;
    T* auxiliar_matrix = nullptr;
    T value = min;
    T* buffer = nullptr;
    T* denominator = nullptr;
    T* temp1 = nullptr;
    T* temp2 = nullptr;
    
    try {
        diagonal = new T[order];
        auxiliar_matrix = new T[order * order];
        if (order > 64)
            {
            buffer = new T[order * order];
            denominator = new T[order * order];
            temp1 = new T[order];
            temp2 = new T[order];
            }
        }
    catch (int)
        {
        if (diagonal != nullptr)
            delete[] diagonal;
        if (auxiliar_matrix != nullptr)
            delete[] auxiliar_matrix;
        if (buffer != nullptr)
            delete[] buffer;
        if (denominator != nullptr)
            delete[] denominator;
        if (temp1 != nullptr)
            delete[] temp1;
        if (temp2 != nullptr)
            delete[] temp2;
        determinant_exception_handle = 1;
        return(-1);
        }
    for (i = 0; i < (order * order); i++) // Copying the matrix
        auxiliar_matrix[i] = pointer[i];
    
    for (i = 0; i < order; i++) // Copying the diagonal of the matrix
        diagonal[i] = auxiliar_matrix[i * (order + 1)];
    
    for (i = 0; i < count; i++) // Calculating the determinant of te values array
        {
        for (j = 0; j < order; j++)
            auxiliar_matrix[j * (order + 1)] =  value + diagonal[j];
            
        output_values[i] = this->Determinant(order, auxiliar_matrix, buffer, denominator, temp1, temp2);
        value = value + step;
        }
    delete[] diagonal;
    delete[] auxiliar_matrix;
    
    if (order > 64)
        {
        delete[] buffer;
        delete[] denominator;
        delete[] temp1;
        delete[] temp2;
        }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::basis_set_Determinant_solver(unsigned int order, T* pointer)
    {   // multithreading code
    unsigned int i, j;
    unsigned int count_eigenvalues = 0;
    T* determinant_values = nullptr;
    T row, max_row, range_coefficient;
    T max_threshold = 0.01;
        
    max_row = 0;
    for (i = 0; i < order; i++)
        {
        row = 0;
        for (j = 0; j < order; j++)
            {
            if (j != i)
                row = row + pointer[(i * order) + j];
            }
        if (abs(row) > max_row)
            max_row = abs(row);
        }
    range_coefficient = (max_row/order) * 2;
    try {
        determinant_values = new T[10000 * order];
        } // Calculating of determinant values array
    catch (int){
        determinant_exception_handle = 1;
        return (-1);
        }
    # pragma omp parallel
        {
        # pragma omp for
        for (i = 0; i < order * 4; i++)
            {
            basis_set_Determinant_set(order, pointer, 1250, range_coefficient * (0.25 * i - order),
            range_coefficient * 0.0002, determinant_values + 1250 * i);
            }
        }
    // end of multithreading code
    determinants.clear();
    while (count_eigenvalues == 0 and max_threshold < 1)
        {
        for (i =1; i < 10000 * order; i++)
            {
            if ((determinant_values[i] * determinant_values[i] <= max_threshold) and (determinant_values[i] * determinant_values[i] <
            determinant_values[i + 1] * determinant_values[i + 1]) and (determinant_values[i] * determinant_values[i] <
            determinant_values[i - 1] * determinant_values[i - 1]))
                determinants.push_back(range_coefficient * ((double(i)/5000) -1 * order));
            }
        max_threshold = max_threshold * 2;
        count_eigenvalues = determinants.size();
        }
    delete[] determinant_values;
    return(0);
    }
// End of Section 1 - work fith Fock matrix - inherited from Huckel_calculations
// Section 2 - generating the wavefunctions
template <typename T>
int Slater_basis_set_calculations<T>::Wavefunction_lenghts_generate(T* lenghts, unsigned int lenght_order)
    { // Generating lenghts 3D cubes for wavefunction calculations
    unsigned int side, x, y, z, index;
    int x_i, y_i, z_i;
    T distance;
    
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = sqrt((x_i * x_i) + (y_i * y_i) + (z_i * z_i));
                lenghts[index] = distance;
                }
            }
        }
    return(0);
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_1s_generate(T* wavefunction, T* lenghts, int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 1s orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/sqrt(Pi) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * exp(-ro/2);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_2s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 2s orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(2 * sqrt(2 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (2.00 - ro) * exp(-ro/2);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_3s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 3s orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(9 * sqrt(3 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (6.00 - (6 * ro) + (ro * ro)) * exp(-ro/2);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4s orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(192 * sqrt(Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (24.00 - (36 * ro) + 12 * (ro * ro)
                    - (ro * ro * ro)) * exp(-ro/2);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5s orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(600 * sqrt(5 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (120 - (240 * ro) + 120 * (ro * ro) 
                    - 20 * (ro * ro * ro) + (ro * ro * ro * ro)) * exp(-ro/2);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6s orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(4320 * sqrt(6 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (720 - (1800 * ro) + 1200 * (ro * ro) 
                    - 300 * (ro * ro * ro) + 30 * (ro * ro * ro * ro) - (ro * ro * ro * ro * ro)) * exp(-ro/2);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 7s orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(15280 * sqrt(7 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (5040 - (15120 * ro) + 12600 * (ro * ro) 
                    - 4200 * (ro * ro * ro) + 630 * (ro * ro * ro * ro) - 42 * 
                    (ro * ro * ro * ro * ro) + (ro * ro * ro * ro * ro * ro)) * exp(-ro/2);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_2px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 2px orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(2 * sqrt(2 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*cos(theta)
                    wavefunction[index] = const_part * ro * exp(-ro/2) * x_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_2pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 2pz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(2 * sqrt(2 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta)
                    wavefunction[index] = const_part * ro * exp(-ro/2) * z_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_2py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
    unsigned int lenght_order)
    { // Generating wavefunction array for 2py orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(2 * sqrt(2 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*sin(theta)
                    wavefunction[index] = const_part * ro * exp(-ro/2) * y_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_3px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 3px orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(9 * sqrt(Pi/2)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*cos(theta)
                    wavefunction[index] = const_part * (4 * ro - (ro * ro)) * exp(-ro/2) * x_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_3pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order) 
    { // Generating wavefunction array for 3pz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(9 * sqrt(Pi/2)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta)
                    wavefunction[index] = const_part * (4 * ro - (ro * ro)) * exp(-ro/2) * z_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_3py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 3py orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(9 * sqrt(Pi/2)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*sin(theta)
                    wavefunction[index] = const_part * (4 * ro - (ro * ro)) * exp(-ro/2) * y_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4px orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(64 * sqrt(5 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*cos(theta)
                    wavefunction[index] = const_part * ro * (20 - (10 * ro) + (ro * ro)) * exp(-ro/2) * x_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4pz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(64 * sqrt(5 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2 * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta)
                    wavefunction[index] = const_part * ro * (20 - (10 * ro) + (ro * ro)) * exp(-ro/2) * z_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
            else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4py orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(64 * sqrt(5 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*sin(theta)
                    wavefunction[index] = const_part * ro * (20 - (10 * ro) + (ro * ro)) * exp(-ro/2) * y_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5px orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(300 * sqrt(10 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*cos(theta)
                    wavefunction[index] = const_part * ro * (120 - (90 * ro) + (18 * ro * ro)
                    - (ro * ro * ro)) * exp(-ro/2) * x_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5pz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(300 * sqrt(10 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta)
                    wavefunction[index] = const_part * ro * (120 - (90 * ro) + (18 * ro * ro)
                    - (ro * ro * ro)) * exp(-ro/2) * z_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5py orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(300 * sqrt(10 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*sin(theta)
                    wavefunction[index] = const_part * ro * (120 - (90 * ro) + (18 * ro * ro)
                    - (ro * ro * ro)) * exp(-ro/2) * y_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6px orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(864 * sqrt(70 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*cos(theta)
                    wavefunction[index] = const_part * ro * (840 - (840 * ro) + (250 * ro * ro) 
                    - (28 * ro * ro * ro) + (ro * ro * ro * ro)) * exp(-ro/2) * x_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6pz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(864 * sqrt(70 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta)
                    wavefunction[index] = const_part * ro * (840 - (840 * ro) + (250 * ro * ro) 
                    - (28 * ro * ro * ro) + (ro * ro * ro * ro)) * exp(-ro/2) * z_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6py orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(864 * sqrt(70 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*sin(theta)
                    wavefunction[index] = const_part * ro * (840 - (840 * ro) + (250 * ro * ro) 
                    - (28 * ro * ro * ro) + (ro * ro * ro * ro)) * exp(-ro/2) * y_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 7px orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(23520 * sqrt(7 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*cos(theta)
                    wavefunction[index] = const_part * ro * (6720 - (8400 * ro) 
                    + (3360 * ro * ro) - (420 * ro * ro * ro) + (40 * ro * ro * ro * ro)
                    - (ro * ro * ro * ro * ro)) * exp(-ro/2) * x_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 7pz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(23520 * sqrt(7 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta)
                    wavefunction[index] = const_part * ro * (6720 - (8400 * ro) 
                    + (3360 * ro * ro) - (420 * ro * ro * ro) + (40 * ro * ro * ro * ro)
                    - (ro * ro * ro * ro * ro)) * exp(-ro/2) * z_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 7py orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(23520 * sqrt(7 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin(fi)*sin(theta)
                    wavefunction[index] = const_part * ro * (6720 - (8400 * ro) 
                    + (3360 * ro * ro) - (420 * ro * ro * ro) + (40 * ro * ro * ro * ro)
                    - (ro * ro * ro * ro * ro)) * exp(-ro/2) * y_i/distance;
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_3dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 3d(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(36 * sqrt(2 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin^2(theta) * cos(2 fi)
                    wavefunction[index] = const_part * (ro * ro * exp(-ro/2))  * ((x_i * x_i) - (y_i * y_i))/(distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_3dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 3dxz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(36 * sqrt(Pi/2)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta) * sin(theta) * cos(fi)
                    wavefunction[index] = const_part * (ro * ro * exp(-ro/2) * (z_i/distance) * x_i/distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_3dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 3dz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(36 * sqrt(6 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (ro * ro * exp(-ro/2)) * 
                    ((3 * z_i * z_i) - (distance * distance))/(distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_3dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order) 
    {// Generating wavefunction array for 3dyz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(36 * sqrt(Pi/2)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta) * sin(theta) * sin(fi)
                    wavefunction[index] = const_part * (ro * ro * exp(-ro/2) * (z_i/distance) * y_i/distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_3dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 3dxy orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(36 * sqrt(2 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin^2(theta) * sin(2 fi)
                    wavefunction[index] = const_part * (ro * ro * exp(-ro/2) * (2 * (y_i/distance)) * (x_i/distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4d(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(384 * sqrt(Pi/3)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin^2(theta) * cos(2 fi)
                    wavefunction[index] = const_part * (6 - ro) * (ro * ro * exp(-ro/2)) 
                    * ((x_i * x_i) - (y_i * y_i))/(distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4dxz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(384 * sqrt(Pi/12)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta) * sin(theta) * cos(fi)
                    wavefunction[index] = const_part * (6 - ro) * (ro * ro * exp(-ro/2) * (z_i/distance) * x_i/distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4dz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(384 * sqrt(Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (6 - ro) * (ro * ro * exp(-ro/2)) * 
                    ((3 * z_i * z_i) - (distance * distance))/(distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order) 
    { // Generating wavefunction array for 4dyz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(384 * sqrt(Pi/12)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta) * sin(theta) * sin(fi)
                    wavefunction[index] = const_part * (6 - ro) * (ro * ro * exp(-ro/2) 
                    * (z_i/distance) * y_i/distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4dxy orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(384 * sqrt(Pi/3)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin^2(theta) * sin(2 fi)
                    wavefunction[index] = const_part * (6 - ro) * (ro * ro * exp(-ro/2) * (2 * (y_i/distance)) * (x_i/distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5d(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(600 * sqrt(70/15 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin^2(theta) * cos(2 fi)
                    wavefunction[index] = const_part * (42 - (14 * ro) + (ro * ro)) * (ro * ro * exp(-ro/2)) 
                    * ((x_i * x_i) - (y_i * y_i))/(distance * distance);
                     norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                   wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5dxz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(600 * sqrt(7/6 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta) * sin(theta) * cos(fi)
                    wavefunction[index] = const_part * (42 - (14 * ro) + (ro * ro)) * (ro * ro * exp(-ro/2) * 
                    (z_i/distance) * x_i/distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5dz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(600 * sqrt(14 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (42 - (14 * ro) + (ro * ro)) * (ro * ro * exp(-ro/2)) * 
                    ((3 * z_i * z_i) - (distance * distance))/(distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order) 
    { // Generating wavefunction array for 5dyz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(600 * sqrt(7/6 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta) * sin(theta) * sin(fi)
                    wavefunction[index] = const_part * (42 - (14 * ro) + (ro * ro)) * (ro * ro * exp(-ro/2) 
                    * (z_i/distance) * y_i/distance);
                     norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5dxy orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(600 * sqrt(70/15 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin^2(theta) * sin(2 fi)
                    wavefunction[index] = const_part * (42 - (14 * ro) + (ro * ro)) 
                    * (ro * ro * exp(-ro/2) * (2 * (y_i/distance)) * (x_i/distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6d(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(3456 * sqrt(7 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin^2(theta) * cos(2 fi)
                    wavefunction[index] = const_part * (336 - (168 * ro) + (24 * ro * ro) - (ro * ro * ro)) 
                    * (ro * ro * exp(-ro/2)) * ((x_i * x_i) - (y_i * y_i))/(distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6dxz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(3456 * sqrt(1.75 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta) * sin(theta) * cos(fi)
                    wavefunction[index] = const_part * (336 - (168 * ro) + (24 * ro * ro) - (ro * ro * ro)) 
                    * (ro * ro * exp(-ro/2) * (z_i/distance) * x_i/distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6dz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(3456 * sqrt(21 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (336 - (168 * ro) + (24 * ro * ro) - (ro * ro * ro)) 
                    * (ro * ro * exp(-ro/2)) * ((3 * z_i * z_i) - (distance * distance))/(distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order) 
    { // Generating wavefunction array for 6dyz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(3456 * sqrt(1.75 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta) * sin(theta) * sin(fi)
                    wavefunction[index] = const_part * (336 - (168 * ro) + (24 * ro * ro) - (ro * ro * ro)) 
                    * (ro * ro * exp(-ro/2) * (z_i/distance) * y_i/distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6dxy orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(3456 * sqrt(7 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin^2(theta) * sin(2 fi)
                    wavefunction[index] = const_part * (336 - (168 * ro) + (24 * ro * ro) - (ro * ro * ro)) 
                    * (ro * ro * exp(-ro/2) * (2 * (y_i/distance)) * (x_i/distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6d(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(28224 * sqrt(7 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin^2(theta) * cos(2 fi)
                    wavefunction[index] = const_part * (3024 - (2016 * ro) + (432 * ro * ro) - (36 * ro * ro * ro) + (ro * ro * ro * ro)) 
                    * (ro * ro * exp(-ro/2)) * ((x_i * x_i) - (y_i * y_i))/(distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6dxz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(28224 * sqrt(1.75 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta) * sin(theta) * cos(fi)
                    wavefunction[index] = const_part * (3024 - (2016 * ro) + (432 * ro * ro) - (36 * ro * ro * ro) + (ro * ro * ro * ro)) 
                    * (ro * ro * exp(-ro/2) * (z_i/distance) * x_i/distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6dz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(28224 * sqrt(21 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (3024 - (2016 * ro) + (432 * ro * ro) - (36 * ro * ro * ro) + (ro * ro * ro * ro)) 
                    * (ro * ro * exp(-ro/2)) * ((3 * z_i * z_i) - (distance * distance))/(distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order) 
    { // Generating wavefunction array for 6dyz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(28224 * sqrt(1.75 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // cos(theta) * sin(theta) * sin(fi)
                    wavefunction[index] = const_part * (3024 - (2016 * ro) + (432 * ro * ro) - (36 * ro * ro * ro) + (ro * ro * ro * ro)) 
                    * (ro * ro * exp(-ro/2) * (z_i/distance) * y_i/distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 6dxy orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(28224 * sqrt(7 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // sin^2(theta) * sin(2 fi)
                    wavefunction[index] = const_part * (3024 - (2016 * ro) + (432 * ro * ro) - (36 * ro * ro * ro) + (ro * ro * ro * ro))
                    * (ro * ro * exp(-ro/2) * (2 * (y_i/distance)) * (x_i/distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4fx(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(384 * sqrt(2 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (ro * ro * ro * exp(-ro/2)) * 
                    (x_i * (x_i * x_i - 3 * y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4fz(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(384 * sqrt(Pi/3)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (ro * ro * ro * exp(-ro/2)) * 
                    (z_i * (x_i * x_i - y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4fxz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(768 * sqrt(Pi/1.2)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (ro * ro * ro * exp(-ro/2)) * 
                    (x_i * ((5 * z_i * z_i) - (distance * distance))/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4fz3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(384 * sqrt(Pi/0.2)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * exp(-ro/2)) * 
                    (z_i * ((5 * z_i * z_i) - (3 * distance * distance))/
                    (distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4fxz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(768 * sqrt(Pi/1.2)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (ro * ro * ro * exp(-ro/2)) * 
                    (y_i * ((5 * z_i * z_i) - (distance * distance))/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4fxyz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(384 * sqrt(Pi/3)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * exp(-ro/2)) * 
                    (2 * x_i * y_i * z_i)/(distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];} // 3 cos(theta)^2 -1
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_4fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 4fx(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(384 * sqrt(2 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/2  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (ro * ro * ro * exp(-ro/2)) * 
                    (y * (3 * x_i * x_i - y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fx(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(2400 * sqrt(Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (8 - ro) * (ro * ro * ro * exp(-ro/2)) * 
                    (x_i * (x_i * x_i - 3 * y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fz(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(1200 * sqrt(Pi/1.5)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (8 - ro) * (ro * ro * ro * exp(-ro/2)) * 
                    (z_i * (x_i * x_i - y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fxz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(2400 * sqrt(Pi/0.6)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5 * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (8 - ro) * (ro * ro * ro * exp(-ro/2)) * 
                    (x_i * ((5 * z_i * z_i) - (distance * distance))/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fz3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(1200 * sqrt(Pi/0.1)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (8 - ro) * (ro * ro * ro * exp(-ro/2)) * 
                    (z_i * ((5 * z_i * z_i) - (3 * distance * distance))/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fxz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(2400 * sqrt(Pi/0.6)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (8 - ro) * (ro * ro * ro * exp(-ro/2)) * 
                    (y_i * ((5 * z_i * z_i) - (distance * distance))/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fxyz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(1200 * sqrt(Pi/1.5)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (8 - ro) * (ro * ro * ro * exp(-ro/2)) * 
                    (2 * x_i * y_i * z_i)/(distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fx(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(2400 * sqrt(Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (8 - ro) * (ro * ro * ro * exp(-ro/2)) * 
                    (y_i * (3 * x_i * x_i - y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fx(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(20736 * sqrt(Pi/2.0)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/6  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (72 - (18 * ro) + (ro * ro)) * (ro * ro * ro * exp(-ro/2)) * 
                    (x_i * (x_i * x_i - 3 * y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fz(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(10368 * sqrt(Pi/3.0)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/6  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (72 - (18 * ro) + (ro * ro)) * (ro * ro * ro * exp(-ro/2)) * 
                    (z_i * (x_i * x_i - y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fxz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(20736 * sqrt(Pi/1.2)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/6  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (72 - (18 * ro) + (ro * ro)) * (ro * ro * ro * exp(-ro/2)) * 
                    (x_i * ((5 * z_i * z_i) - (distance * distance))/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fz3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(10368 * sqrt(Pi/0.2)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/6  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (72 - (18 * ro) + (ro * ro)) * (ro * ro * ro * exp(-ro/2)) * 
                    (z_i * ((5 * z_i * z_i) - (3 * distance * distance))/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fxz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(20736 * sqrt(Pi/1.2)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/6  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (72 - (18 * ro) + (ro * ro)) * (ro * ro * ro * exp(-ro/2)) * 
                    (y_i * ((5 * z_i * z_i) - (distance * distance))/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fxyz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(10368 * sqrt(Pi/3.0)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/6  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (72 - (18 * ro) + (ro * ro)) * (ro * ro * ro * exp(-ro/2)) * 
                    (2 * x_i * y_i * z_i)/(distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fx(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(20736 * sqrt(Pi/2.0)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/6  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (72 - (18 * ro) + (ro * ro)) * (ro * ro * ro * exp(-ro/2)) * 
                    (y_i * (3 * x_i * x_i - y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fx(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(141120 * sqrt(Pi/2.4)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (720 - (270 * ro) + (30 * ro * ro) - (ro * ro * ro))
                    * (ro * ro * ro * exp(-ro/2)) * (x_i * (x_i * x_i - 3 * y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fz(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(70560 * sqrt(Pi/3.6)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (720 - (270 * ro) + (30 * ro * ro) - (ro * ro * ro))
                    * (ro * ro * ro * exp(-ro/2)) * (z_i * (x_i * x_i - y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fxz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(141120 * sqrt(Pi/1.44)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (720 - (270 * ro) + (30 * ro * ro) - (ro * ro * ro))
                    * (ro * ro * ro * exp(-ro/2)) * (x_i * ((5 * z_i * z_i) - (distance * distance))/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fz3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(70560 * sqrt(Pi/0.24)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (720 - (270 * ro) + (30 * ro * ro) - (ro * ro * ro))* (ro * ro * ro * exp(-ro/2))
                    * (z_i * ((5 * z_i * z_i) - (3 * distance * distance))/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fxz2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 1/(141120 * sqrt(Pi/1.44)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (720 - (270 * ro) + (30 * ro * ro) - (ro * ro * ro))
                    * (ro * ro * ro * exp(-ro/2)) * (y_i * ((5 * z_i * z_i) - (distance * distance))/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fxyz orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(70560 * sqrt(Pi/3.6)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (720 - (270 * ro) + (30 * ro * ro) - (ro * ro * ro))
                    * (ro * ro * ro * exp(-ro/2)) * (2 * x_i * y_i * z_i)/(distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
unsigned int lenght_order)
    { // Generating wavefunction array for 5fx(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(141120 * sqrt(Pi/2.4)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){ // 3 cos(theta)^2 -1
                    wavefunction[index] = const_part * (720 - (270 * ro) + (30 * ro * ro) - (ro * ro * ro))
                    * (ro * ro * ro * exp(-ro/2)) * (y_i * (3 * x_i * x_i - y_i * y_i)/(distance * distance * distance));
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5gz4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 5gz4 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(900 * sqrt(630/256 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * exp(-ro/2)) * (35 * z_i * z_i * z_i * z_i
                    - 30 * z_i * z_i * distance * distance + 3 * distance * distance * distance)
                    /(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5gz3y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 5gz3y orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(900 * sqrt(98.4375 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * exp(-ro/2)) * y_i * z_i *
                    (7 * z_i * z_i - 3 * distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5gz3x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 5gz3x orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(900 * sqrt(98.4375 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * exp(-ro/2)) * x_i * z_i *
                    (7 * z_i * z_i - 3 * distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5gz2xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 5gz2xy orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(900 * sqrt(49.21875 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * exp(-ro/2)) * 2 * x_i * y_i *
                    (7 * z_i * z_i - distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5gz2_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 5gz2(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(900 * sqrt(49.21875 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * exp(-ro/2)) * (x_i * x_i - y_i * y_i) *
                    (7 * z_i * z_i - distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5gz_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 5gzy3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(900 * sqrt(689.0625 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * exp(-ro/2)) * y_i * z_i *
                    (3 * x_i * x_i - y_i * y_i)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5gz_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 5gzx3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(900 * sqrt(689.0625 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * exp(-ro/2)) * x_i * z_i *
                    (x_i * x_i - 3 * y_i * y_i)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5gxy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 5gxy(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(900 * sqrt(86.1328125 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * exp(-ro/2)) * 4 * x_i * y_i *
                    (x_i * x_i - y_i * y_i)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_5gx4_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 5gxy(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(900 * sqrt(86.1328125 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/5  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * exp(-ro/2)) *
                    (x_i * x_i * x_i * x_i + y_i * y_i * y_i * y_i - 6 * x_i * x_i * y_i * y_i)/
                    (distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6gz4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6gz4 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(12960 * sqrt(63/256 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (10 - ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    (35 * z_i * z_i * z_i * z_i - 30 * z_i * z_i * distance * distance +
                    3 * distance * distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6gz3y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6gz3y orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(12960 * sqrt(9.84375 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (10 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * y_i * z_i *
                    (7 * z_i * z_i - 3 * distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6gz3x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6gz3x orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(12960 * sqrt(9.84375 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (10 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * x_i * z_i *
                    (7 * z_i * z_i - 3 * distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6gz2xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6gz2xy orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(12960 * sqrt(4.921875 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (10 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * 2 * x_i * y_i *
                    (7 * z_i * z_i - distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6gz2_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6gz2(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(12960 * sqrt(4.921875 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (10 - ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    (x_i * x_i - y_i * y_i) * (7 * z_i * z_i - distance * distance)/
                    (distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6gz_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6gzy3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(12960 * sqrt(68.90625 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (10 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * y_i * z_i *
                    (3 * x_i * x_i - y_i * y_i)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6gz_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6gzx3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(12960 * sqrt(68.90625 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (10 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * x_i * z_i *
                    (x_i * x_i - 3 * y_i * y_i)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6gxy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6gxy(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(12960 * sqrt(8.61328125 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (10 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * 4 * x_i * y_i *
                    (x_i * x_i - y_i * y_i)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6gx4_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6gxy(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(12960 * sqrt(8.61328125 * Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (10 - ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    (x_i * x_i * x_i * x_i + y_i * y_i * y_i * y_i - 6 * x_i * x_i * y_i * y_i)/
                    (distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7gz4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7gz4 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(17640 * sqrt(13860/256 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (110 - 22 * ro + ro * ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    (35 * z_i * z_i * z_i * z_i - 30 * z_i * z_i * distance * distance +
                    3 * distance * distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7gz3y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7gz3y orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(17640 * sqrt(216.5625 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (110 - 22 * ro + ro * ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    y_i * z_i * (7 * z_i * z_i - 3 * distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7gz3x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7gz3x orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(17640 * sqrt(216.5625 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (110 - 22 * ro + ro * ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    x_i * z_i * (7 * z_i * z_i - 3 * distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7gz2xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7gz2xy orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(17640 * sqrt(108.28125 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (110 - 22 * ro + ro * ro) * (ro * ro * ro * ro * exp(-ro/2)) * 2 *
                    x_i * y_i * (7 * z_i * z_i - distance * distance)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7gz2_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7gz2(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(17640 * sqrt(108.28125 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (110 - 22 * ro + ro * ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    (x_i * x_i - y_i * y_i) * (7 * z_i * z_i - distance * distance)/
                    (distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7gz_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7gzy3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(17640 * sqrt(1515.9375 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (110 - 22 * ro + ro * ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    y_i * z_i * (3 * x_i * x_i - y_i * y_i)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7gz_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7gzx3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(17640 * sqrt(1515.9375 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (110 - 22 * ro + ro * ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    x_i * z_i * (x_i * x_i - 3 * y_i * y_i)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7gxy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7gxy(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(17640 * sqrt(189.4921875 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (110 - 22 * ro + ro * ro) * (ro * ro * ro * ro * exp(-ro/2)) * 4 *
                    x_i * y_i * (x_i * x_i - y_i * y_i)/(distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7gx4_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7gxy(x2-y2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part =1/(17640 * sqrt(189.4921875 *Pi)) * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2*Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (110 - 22 * ro + ro * ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    (x_i * x_i * x_i * x_i + y_i * y_i * y_i * y_i - 6 * x_i * x_i * y_i * y_i)/
                    (distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6hz5_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6hz5 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 5.87591177054549E-09 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part  * (ro * ro * ro * ro * ro * exp(-ro/2)) *
                    z_i * (63 * z_i * z_i * z_i * z_i - 70 * z_i * z_i * distance * distance +
                    15 * distance * distance * distance * distance)/
                    (distance * distance * distance * distance  * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6hz4y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6hz4y orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 5.87591177054549E-09 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * exp(-ro/2)) * y_i *
                    (21 * z_i * z_i * z_i * z_i - 14 * z_i * z_i * distance * distance +
                    distance * distance * distance * distance)/(distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6hz4x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6hz4x orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 5.87591177054549E-09 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * exp(-ro/2)) * x_i *
                    (21 * z_i * z_i * z_i * z_i - 14 * z_i * z_i * distance * distance +
                    distance * distance * distance * distance)/(distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6hz3xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6hz3xy orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 5.87591177054549E-09 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * exp(-ro/2)) *
                    (2 * x_i * y_i) * (3 * z_i * z_i * z_i - z_i * distance * distance)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6hz3_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6hz3_x2-y2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 5.87591177054549E-09 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * exp(-ro/2)) *
                    (x_i * x_i - y_i * y_i) * (3 * z_i * z_i * z_i - distance * z_i * z_i)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6hz2_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6hz2y3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 5.87591177054549E-09 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * exp(-ro/2)) * y_i * 
                    ( 3 * x_i * x_i - y_i * y_i) * (9 * z_i * z_i - distance * distance)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6hz2_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6hz2_x3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 5.87591177054549E-09 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * exp(-ro/2)) * x_i * 
                    (x_i * x_i - 3 * y_i * y_i) * (9 * z_i * z_i - distance * distance)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6hz_x3y_xy3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6hz(x4-x2y2+y4) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 5.87591177054549E-09 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * exp(-ro/2)) * z_i * 
                    (4 * x_i * x_i * x_i * y_i - 4 * x_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6hz_x4_x2y2_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6hz(x4-x2y2+y4) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 5.87591177054549E-09 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * exp(-ro/2)) * z_i * 
                    (x_i * x_i * x_i * x_i - 6 * x_i * x_i * y_i * y_i + y_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6hy_x4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6hyx4 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 5.87591177054549E-09 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * exp(-ro/2)) * y_i * 
                    (5 * x_i * x_i * x_i * x_i - 10 * x_i * x_i * y_i * y_i + y_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_6hx_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 6hy4x orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 5.87591177054549E-09 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = Z/3  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * exp(-ro/2)) * x_i * 
                    (x_i * x_i * x_i * x_i - 10 * x_i * x_i * y_i * y_i + 5 * y_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7hz5_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hz5 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 7.26569729185807E-10 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (12 - ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    z_i * (63 * z_i * z_i * z_i * z_i - 70 * z_i * z_i * distance * distance +
                    15 * distance * distance * distance * distance)/
                    (distance * distance * distance * distance  * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7hz4y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hz4y orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 7.26569729185807E-10 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7 * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (12 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * y_i *
                    (21 * z_i * z_i * z_i * z_i - 14 * z_i * z_i * distance * distance +
                    distance * distance * distance * distance)/(distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7hz4x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hz4x orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 7.26569729185807E-10 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (12 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * x_i *
                    (21 * z_i * z_i * z_i * z_i - 14 * z_i * z_i * distance * distance +
                    distance * distance * distance * distance)/(distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7hz3xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hz3xy orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 7.26569729185807E-10 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (12 - ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    (2 * x_i * y_i) * (3 * z_i * z_i * z_i - z_i * distance * distance)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7hz3_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hz3_x2-y2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 7.26569729185807E-10 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (12 - ro) * (ro * ro * ro * ro * exp(-ro/2)) *
                    (x_i * x_i - y_i * y_i) * (3 * z_i * z_i * z_i - distance * z_i * z_i)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7hz2_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hz2y3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 7.26569729185807E-10 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (12 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * y_i * 
                    ( 3 * x_i * x_i - y_i * y_i) * (9 * z_i * z_i - distance * distance)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7hz2_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hz2_x3 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 7.26569729185807E-10 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (12 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * x_i * 
                    (x_i * x_i - 3 * y_i * y_i) * (9 * z_i * z_i - distance * distance)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7hz_x3y_xy3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hz(x4-x2y2+y4) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 7.26569729185807E-10 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (12 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * z_i * 
                    (4 * x_i * x_i * x_i * y_i - 4 * x_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7hz_x4_x2y2_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hz(x4-x2y2+y4) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 7.26569729185807E-10 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (12 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * z_i * 
                    (x_i * x_i * x_i * x_i - 6 * x_i * x_i * y_i * y_i + y_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7hy_x4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hyx4 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 7.26569729185807E-10 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (12 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * y_i * 
                    (5 * x_i * x_i * x_i * x_i - 10 * x_i * x_i * y_i * y_i + y_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7hx_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hxy4 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 7.26569729185807E-10 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (12 - ro) * (ro * ro * ro * ro * exp(-ro/2)) * x_i * 
                    (x_i * x_i * x_i * x_i - 10 * x_i * x_i * y_i * y_i + 5 * y_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7iz6_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7iz6 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) *
                    (231 * z_i * z_i * z_i * z_i * z_i * z_i - 315 * z_i * z_i * z_i * z_i * distance * distance +
                    105 * z_i * z_i * distance * distance * distance * distance)/
                    (distance * distance * distance * distance  * distance  * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7iz5y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7iz5y orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7 * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) * y_i *
                    (33 * z_i * z_i * z_i * z_i * z_i - 30 * z_i * z_i * z_i * distance * distance +
                    z_i * distance * distance * distance * distance)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7iz5x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7iz5x orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) * x_i *
                    (33 * z_i * z_i * z_i * z_i * z_i - 30 * z_i * z_i * z_i * distance * distance +
                    z_i * distance * distance * distance * distance)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7iz4xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7iz4xy orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) *
                    (2 * x_i * y_i) * (33 * z_i * z_i * z_i * z_i - 18 * z_i * z_i * distance * distance +
                    distance * distance * distance * distance)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7iz4_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7hz4_x2-y2 orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) *
                    (x_i * x_i - y_i * y_i) * (33 * z_i * z_i * z_i * z_i - 18 * z_i * z_i * distance * distance +
                    distance * distance * distance * distance)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7iy__x2_y2__z3_zr2_generate(T* wavefunction, T* lenghts, unsigned int Z,
T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7iy(x2-y2)(z3-zr2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) * y_i * 
                    (x_i * x_i - 3 * y_i * y_i) * (11 * z_i * z_i * z_i - 3 * z_i * distance * distance)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7ix__x2_y2__z3_zr2_generate(T* wavefunction, T* lenghts, unsigned int Z,
T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7ix(x2-y2)(z3-zr2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) * x_i * 
                    (x_i * x_i - 3 * y_i * y_i) * (11 * z_i * z_i * z_i - 3 * z_i * distance * distance)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7i_xy__x2_y2__z2_r2_generate(T* wavefunction, T* lenghts, unsigned int Z,
T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7ii(xy)(x2-y2)(z2-r2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) * x_i * y_i *
                    (4 * x_i * x_i * x_i * x_i - 6 * x_i * x_i * y_i * y_i + y_i * y_i * y_i * y_i) *
                    (11 * z_i * z_i - distance * distance)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7i_x4_x2y2_y4__z2_r2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7i(x4-x2y2+y4)(z2-r2) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) *
                    (x_i * x_i * x_i * x_i - 6 * x_i * x_i * y_i * y_i + y_i * y_i * y_i * y_i) *
                    (11 * z_i * z_i - distance * distance)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7i_yz__x4_x2y2_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7iyz(x4-x2y2+y4) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) * y_i * z_i *
                    (5 * x_i * x_i * x_i * x_i - 10 * x_i * x_i * y_i * y_i + y_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7i_xz__x4_x2y2_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7ixz(x4-x2y2+y4) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        z_i = (int(z) - int(lenght_order));
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) * x_i * z_i *
                    (x_i * x_i * x_i * x_i - 10 * x_i * x_i * y_i * y_i + 5 * y_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7i_yx5_x3y3_xy5_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7i(yx5-x3y3+xy5) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) *
                    (6 * y_i * x_i * x_i * x_i * x_i * x_i - 20 * x_i * x_i * x_i * y_i * y_i * y_i +
                    6 * x_i * y_i * y_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
T Slater_basis_set_calculations<T>::Wavefunction_7i_x6__x4y2_x2y4_y6_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
    { // Generating wavefunction array for 7i(yx5-x3y3+xy5) orbital
    unsigned int side;
    unsigned int index;
    unsigned int x, y, z;
    int x_i, y_i, z_i;
    T ro_0, ro;
    T distance, const_part;
    T norm_coeff = 0;
    
    const_part = 2.87877408107628E-11 * Z * Z * 1/sqrt(Z); // calculating of constant part of wavefunction
    ro_0 = 2 * Z/7  * multiplier; // Algorithm based on Laguerr polynoms
    side = (lenght_order * 2) + 1;
    for (z = 0; z < side; z++)
        {
        for (y = 0; y < side; y++)
            {
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
                index = x + (y * side) + (z * side * side);
                distance = lenghts[index];
                ro = ro_0 * distance;
                if (distance != 0){
                    wavefunction[index] = const_part * (ro * ro * ro * ro * ro * ro * exp(-ro/2)) *
                    (x_i * x_i * x_i * x_i * x_i * x_i - 15 * x_i * x_i * x_i * x_i * y_i * y_i +
                    + 15 * x_i * x_i * y_i * y_i * y_i * y_i - y_i * y_i * y_i * y_i * y_i * y_i)/
                    (distance * distance * distance * distance * distance * distance);
                    norm_coeff = norm_coeff + wavefunction[index] * wavefunction[index];}
                else
                    wavefunction[index] = 0;
                }
            }
        }
    return(sqrt(norm_coeff));
    }
template <typename T>
int Slater_basis_set_calculations<T>::Wavefunction_normalize(T* wavefunction_pointer, T normalisation_constant, unsigned int size)
    {
    unsigned int i;
    
    for (i = 0; i < size - 7; i = i + 8)
        { // Vectorisation code
        wavefunction_pointer[i] = wavefunction_pointer[i]/normalisation_constant; 
        wavefunction_pointer[i + 1] = wavefunction_pointer[i + 1]/normalisation_constant;
        wavefunction_pointer[i + 2] = wavefunction_pointer[i + 2]/normalisation_constant;
        wavefunction_pointer[i + 3] = wavefunction_pointer[i + 3]/normalisation_constant;
        wavefunction_pointer[i + 4] = wavefunction_pointer[i + 4]/normalisation_constant;
        wavefunction_pointer[i + 5] = wavefunction_pointer[i + 5]/normalisation_constant;
        wavefunction_pointer[i + 6] = wavefunction_pointer[i + 6]/normalisation_constant;
        wavefunction_pointer[i + 7] = wavefunction_pointer[i + 7]/normalisation_constant;
        } // End of vectorisation code
    for (i = size - (size % 8); i < size; i++)
        wavefunction_pointer[i] = wavefunction_pointer[i]/normalisation_constant; // Finishing calculation
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Orbitals_to_wavefunctions(unsigned int n, unsigned int l, int m, unsigned int lenght_order, T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int* x_range,
unsigned int* y_range, unsigned int* z_range)
    {
    unsigned int size;
    T normalisation_constant = 0;
    
    if (n == 1)
        normalisation_constant = Wavefunction_1s_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
        
    if (n == 2)
        {
        if (l == 0)
            normalisation_constant = Wavefunction_2s_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            
        if (l == 1)
            {
                if (m == 1)
                    normalisation_constant = Wavefunction_2px_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
                if (m == 0)
                    normalisation_constant = Wavefunction_2pz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
                if (m == -1)
                    normalisation_constant = Wavefunction_2py_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        }
    if (n == 3)
        {
        if (l == 0)
            normalisation_constant = Wavefunction_3s_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            
        if (l == 1)
            {
            if (m == 1)
                normalisation_constant = Wavefunction_3px_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_3pz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_3py_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 2)
            {
            if (m == 2)
                normalisation_constant = Wavefunction_3dx2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_3dxz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_3dz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_3dyz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_3dxy_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        }
    if (n == 4)
        {
        if (l == 0)
            normalisation_constant = Wavefunction_4s_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            
        if (l == 1)
            {
            if (m == 1)
                normalisation_constant = Wavefunction_4px_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_4pz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_4py_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 2)
            {
            if (m == 2)
                normalisation_constant = Wavefunction_4dx2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_4dxz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_4dz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_4dyz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_4dxy_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 3)
            {
            if (m == 3)
                normalisation_constant = Wavefunction_4fx_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 2)
                normalisation_constant = Wavefunction_4fz_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_4fxz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_4fz3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_4fyz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_4fxyz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -3)
                normalisation_constant = Wavefunction_4fy_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        }
    if (n == 5)
        {
        if (l == 0)
            normalisation_constant = Wavefunction_5s_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            
        if (l == 1)
            {
            if (m == 1)
                normalisation_constant = Wavefunction_5px_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_5pz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_5py_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 2)
            {
            if (m == 2)
                normalisation_constant = Wavefunction_5dx2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_5dxz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_5dz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_5dyz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_5dxy_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 3)
            {
            if (m == 3)
                normalisation_constant = Wavefunction_5fx_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 2)
                normalisation_constant = Wavefunction_5fz_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_5fxz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_5fz3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_5fyz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_5fxyz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -3)
                normalisation_constant = Wavefunction_5fy_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 4)
            {
            if (m == 4)
                normalisation_constant = Wavefunction_5gz4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 3)
                normalisation_constant = Wavefunction_5gz3y_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 2)
                normalisation_constant = Wavefunction_5gz3x_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_5gz2xy_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_5gz2_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_5gz_x3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_5gz_x3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -3)
                normalisation_constant = Wavefunction_5gxy_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -4)
                normalisation_constant = Wavefunction_5gx4_y4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        }
    if (n == 6)
        {
        if (l == 0)
            normalisation_constant = Wavefunction_6s_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            
        if (l == 1)
            {
            if (m == 1)
                normalisation_constant = Wavefunction_6px_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_6pz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_6py_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 2)
            {
            if (m == 2)
                normalisation_constant = Wavefunction_6dx2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_6dxz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_6dz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_6dyz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_6dxy_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 3)
            {
            if (m == 3)
                normalisation_constant = Wavefunction_6fx_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 2)
                normalisation_constant = Wavefunction_6fz_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_6fxz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_6fz3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_6fyz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_6fxyz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -3)
                normalisation_constant = Wavefunction_6fy_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 4)
            {
            if (m == 4)
                normalisation_constant = Wavefunction_6gz4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 3)
                normalisation_constant = Wavefunction_6gz3y_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 2)
                normalisation_constant = Wavefunction_6gz3x_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_6gz2xy_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_6gz2_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_6gz_x3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_6gz_x3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -3)
                normalisation_constant = Wavefunction_6gxy_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -4)
                normalisation_constant = Wavefunction_6gx4_y4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 5)
            {
            if (m == 5)
                normalisation_constant = Wavefunction_6hz5_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 4)
                normalisation_constant = Wavefunction_6hz4y_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 3)
                normalisation_constant = Wavefunction_6hz4x_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 2)
                normalisation_constant = Wavefunction_6hz3xy_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_6hz3_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_6hz2_y3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_6hz2_x3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_6hz_x3y_xy3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -3)
                normalisation_constant = Wavefunction_6hz_x4_x2y2_y4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -4)
                normalisation_constant = Wavefunction_6hy_x4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -5)
                normalisation_constant = Wavefunction_6hx_y4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        }
    if (n == 7)
        {
        if (l == 0)
            normalisation_constant = Wavefunction_7s_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            
        if (l == 1)
            {
            if (m == 1)
                normalisation_constant = Wavefunction_7px_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_7pz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_7py_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 2)
            {
            if (m == 2)
                normalisation_constant = Wavefunction_7dx2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_7dxz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_7dz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_7dyz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_7dxy_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 3)
            {
            if (m == 3)
                normalisation_constant = Wavefunction_7fx_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 2)
                normalisation_constant = Wavefunction_7fz_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_7fxz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_7fz3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_7fyz2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_7fxyz_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -3)
                normalisation_constant = Wavefunction_7fy_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 4)
            {
            if (m == 4)
                normalisation_constant = Wavefunction_7gz4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 3)
                normalisation_constant = Wavefunction_7gz3y_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 2)
                normalisation_constant = Wavefunction_7gz3x_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_7gz2xy_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_7gz2_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_7gz_x3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_7gz_x3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -3)
                normalisation_constant = Wavefunction_7gxy_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -4)
                normalisation_constant = Wavefunction_7gx4_y4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 5)
            {
            if (m == 5)
                normalisation_constant = Wavefunction_7hz5_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 4)
                normalisation_constant = Wavefunction_7hz4y_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 3)
                normalisation_constant = Wavefunction_7hz4x_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 2)
                normalisation_constant = Wavefunction_7hz3xy_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_7hz3_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_7hz2_y3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_7hz2_x3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_7hz_x3y_xy3_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -3)
                normalisation_constant = Wavefunction_7hz_x4_x2y2_y4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -4)
                normalisation_constant = Wavefunction_7hy_x4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -5)
                normalisation_constant = Wavefunction_7hx_y4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        if (l == 6)
            {
            if (m == 6)
                normalisation_constant = Wavefunction_7iz6_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 5)
                normalisation_constant = Wavefunction_7iz5y_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 4)
                normalisation_constant = Wavefunction_7iz5x_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 3)
                normalisation_constant = Wavefunction_7iz4xy_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 2)
                normalisation_constant = Wavefunction_7iz4_x2_y2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 1)
                normalisation_constant = Wavefunction_7iy__x2_y2__z3_zr2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == 0)
                normalisation_constant = Wavefunction_7ix__x2_y2__z3_zr2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -1)
                normalisation_constant = Wavefunction_7i_xy__x2_y2__z2_r2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -2)
                normalisation_constant = Wavefunction_7i_x4_x2y2_y4__z2_r2_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -3)
                normalisation_constant = Wavefunction_7i_yz__x4_x2y2_y4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -4)
                normalisation_constant = Wavefunction_7i_xz__x4_x2y2_y4_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -5)
                normalisation_constant = Wavefunction_7i_yx5_x3y3_xy5_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            if (m == -6)
                normalisation_constant = Wavefunction_7i_x6__x4y2_x2y4_y6_generate(wavefunction, lenghts, Z, multiplier, lenght_order);
            }
        }
    size = (2 * lenght_order + 1) * (2 * lenght_order + 1) * (2 * lenght_order + 1);
    Wavefunction_normalize(wavefunction, normalisation_constant, size);
    Wavefunction_range_detect(wavefunction, lenght_order, x_range, y_range, z_range);
    return(0);
    }
template <typename T>
T Slater_basis_set_calculations<T>::Get_relative_Hartree_length(unsigned int Z, unsigned int n)
    { // Setting relative permeability and permeability constant for wavefunction calculations
    T relative_Hartree_lenght;
    T relative_electron_mass;
    T relativistic_shrinkage;
    
    relative_electron_mass = (Z * mp)/(Z * mp + me);
    relativistic_shrinkage = 1/sqrt(1 - pow(Z * hyperfine_structure_constant/n, 2));
    relative_Hartree_lenght = relativistic_shrinkage/relative_electron_mass;
    return(relative_Hartree_lenght);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Wavefunction_range_detect(T* wavefunction_pointer, unsigned int lenght_order,
unsigned int* x_range, unsigned int * y_range, unsigned int * z_range)
    { // Calculate the absolute value of overleap integral for multiplying the energy differences from fock matrix
    unsigned int side;
    unsigned int i, j, k;
    unsigned int i_range, j_range, k_range;
    
    T threshold = wavefunction_integration_threshold;
    
    side = (2 * lenght_order + 1);
    // Is is not a wavefunction integration range detection necessary return max ranges
    if (wavefunction_pointer[0] > threshold)
        {
        x_range[0] = lenght_order;
        y_range[0] = lenght_order;
        z_range[0] = lenght_order;
        return(0);
        }
    
    // detecting a effective x, y and z ranges of integration
    i_range = 0;
    j_range = 0;
    k_range = 0;
    for (k = 0; k < side; k++)
            for (j = 0; j < side; j++)
                for (i = 0; i < side; i++)
                    {
                    if (wavefunction_pointer[i + j * side + k * side * side] > threshold)
                        {
                        if (i > lenght_order)
                            if (i > i_range + lenght_order)
                                i_range = i - lenght_order;
                        else
                            if (lenght_order > i_range + i)
                                i_range = lenght_order - i;
                        if (j > lenght_order)
                            if (j > j_range + lenght_order)
                                j_range = j - lenght_order;
                        else
                            if (lenght_order > j_range + j)
                                j_range = lenght_order - j;
                        if (k > lenght_order)
                            if (k > k_range + lenght_order)
                                k_range = k - lenght_order;
                        else
                            if (lenght_order > k_range + k)
                                k_range = lenght_order - k;
                        }
                    }
    x_range[0] = i_range;
    y_range[0] = j_range;
    z_range[0] = k_range;
    return(0);
    }
// End of Section 2 - generating wavefunctions, Section 3 - mathematical operations for wavefunctions, probabilities densities and
// integrals
template <typename T>
int Slater_basis_set_calculations<T>::Wavefunction_multiply(T* wavefunction_1, T* wavefunction_2, T* probabilities,
unsigned int lenght_order, T d_x, T d_y, T d_z)
    { // Multiplyng of wavefunction
    unsigned int i, j, k;
    unsigned int side;
    unsigned int x_side, y_side, z_side;
    unsigned int x_contraction, y_contraction, z_contraction;
    unsigned int x_1_min, y_1_min, z_1_min;
    unsigned int x_2_min, y_2_min, z_2_min;
    unsigned int x_condition, x_condition_2;
    int x, y, z;
    
    x = d_x * lenght_order/vector_lenght;
    y = d_y * lenght_order/vector_lenght;
    z = d_z * lenght_order/vector_lenght;
    side = (2 * lenght_order + 1);
    x_contraction = abs(x);
    y_contraction = abs(y);
    z_contraction = abs(z);
    
    if (x_contraction > side)
        x_contraction = side;
    if (y_contraction > side)
        y_contraction = side;
    if (z_contraction > side)
        z_contraction = side;
    
    x_side = side - x_contraction;
    y_side = side - y_contraction;
    z_side = side - z_contraction;
    
    if (x_contraction >= side or y_contraction >= side or z_contraction >= side)
        return(-1);
        
    if (x < 0)
        {
        x_1_min = 0;
        x_2_min = x_contraction;
        }
    else
        {
        x_1_min = x_contraction;
        x_2_min = 0;
        }
    if (y < 0)
        {
        y_1_min = 0;
        y_2_min = y_contraction;
        }
    else
        {
        y_1_min = y_contraction;
        y_2_min = 0;
        }
    if (z < 0)
        {
        z_1_min = 0;
        z_2_min = z_contraction;
        }
    else
        {
        z_1_min = z_contraction;
        z_2_min = 0;
        }
    
    x_condition = side - x_contraction;   
    x_condition_2 = x_condition - (x_condition % 8);
    
    // vectorisation code
    for (i = 0; i < (side - z_contraction); i++)
        {
        for (j = 0; j < (side - y_contraction); j++)
            {
            for (k = 0; k + 7 < x_condition; k+=8)
                {
                probabilities[(i * (y_side * x_side)) + (j * x_side) + k] =
                wavefunction_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min)] *
                wavefunction_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min)];
                probabilities[(i * (y_side * x_side)) + (j * x_side) + k + 1] =
                wavefunction_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 1] *
                wavefunction_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 1];
                probabilities[(i * (y_side * x_side)) + (j * x_side) + k + 2] =
                wavefunction_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 2] *
                wavefunction_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 2];
                probabilities[(i * (y_side * x_side)) + (j * x_side) + k + 3] =
                wavefunction_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 3] *
                wavefunction_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 3];
                probabilities[(i * (y_side * x_side)) + (j * x_side) + k + 4] =
                wavefunction_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 4] *
                wavefunction_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 4];
                probabilities[(i * (y_side * x_side)) + (j * x_side) + k + 5] =
                wavefunction_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 5] *
                wavefunction_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 5];
                probabilities[(i * (y_side * x_side)) + (j * x_side) + k + 6] =
                wavefunction_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 6] *
                wavefunction_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 6];
                probabilities[(i * (y_side * x_side)) + (j * x_side) + k + 7] =
                wavefunction_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 7] *
                wavefunction_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 7];
                }
            for (k = x_condition_2; k < x_condition; k++)
                {
                probabilities[(i * (y_side * x_side)) + (j * x_side) + k] =
                wavefunction_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min)] *
                wavefunction_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min)];
                }
            }
        }
        // end of vectorisation code
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Wavefunction_multiply(T* wavefunction_1, T* wavefunction_2, T* probabilities, unsigned int lenght_order)
    {
    // Multiplyng of wavefunctions
    unsigned int side;
    unsigned int size;
    unsigned int i;
    
    side = (2 * lenght_order + 1);
    size = side * side * side;
    // vectorisation code
    for (i = 0; i < (size - 7); i = i + 8) // Multiplying of Wavefunction_1 and Wavefunction_1
        {
        probabilities[i] = wavefunction_1[i] * wavefunction_2[i];
        probabilities[i + 1] = wavefunction_1[i + 1] * wavefunction_2[i + 1];
        probabilities[i + 2] = wavefunction_1[i + 2] * wavefunction_2[i + 2];
        probabilities[i + 3] = wavefunction_1[i + 3] * wavefunction_2[i + 3];
        probabilities[i + 4] = wavefunction_1[i + 4] * wavefunction_2[i + 4];
        probabilities[i + 5] = wavefunction_1[i + 5] * wavefunction_2[i + 5];
        probabilities[i + 6] = wavefunction_1[i + 6] * wavefunction_2[i + 6];
        probabilities[i + 7] = wavefunction_1[i + 7] * wavefunction_2[i + 7]; 
        }
    for (i = size - (size % 8); i < (size); i++)
        probabilities[i] = wavefunction_1[i] * wavefunction_2[i];
    // end of vectorisation code
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Wavefunction_relative_lenghts_generate(T* reverse_relative_lenghts,
unsigned int lenght_order)
    { // Generating lenghts 3D cubes for wavefunction calculations
    unsigned int side, i, j, k, l, m ,n, pre_index, index;
    int x_i, y_i, z_i;
    T distance;
    
    side = (lenght_order * 2) + 1;
    for (i = 0; i < side; i++)
        for (j = 0; j < side; j++)
            for (k = 0; k < side; k++)
                for (l = 0; l < side; l++)
                    {
                    z_i = (int(l) - int(i)) * (int(l) - int(i));
                    pre_index = (i * side * side * side * side * side) +
                    (j * side * side * side * side) + (k * side * side * side) +
                    (l * side * side);
                    for (m = 0; m < side; m++)
                        {
                        y_i = (int(m) - int(j)) * (int(m) - int(j));
                        for (n = 0; n < side; n++)
                            {
                            x_i = (int(n) - int(k)) * (int(n) - int(k));
                            index = pre_index + (m * side) + n;
                            distance = sqrt(x_i + y_i + z_i);
                            if (distance != 0)
                                reverse_relative_lenghts[index] = 1/distance;
                            else
                                reverse_relative_lenghts[index] = 0;
                            }
                        }
                    }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Wavefunction_square(T* wavefunction_1, T* probabilities, unsigned int lenght_order)
    { // Multiplyng of wavefunction
    unsigned int side;
    unsigned int size;
    unsigned int i;
    
    side = (2 * lenght_order + 1);
    size = side * side * side;
    // vectorisation code
    for (i = 0; i < (size - 7); i = i + 8) // Multiplying of Wavefunction_1 and Wavefunction_1
        {
        probabilities[i] = wavefunction_1[i] * wavefunction_1[i];
        probabilities[i + 1] = wavefunction_1[i + 1] * wavefunction_1[i + 1];
        probabilities[i + 2] = wavefunction_1[i + 2] * wavefunction_1[i + 2];
        probabilities[i + 3] = wavefunction_1[i + 3] * wavefunction_1[i + 3];
        probabilities[i + 4] = wavefunction_1[i + 4] * wavefunction_1[i + 4];
        probabilities[i + 5] = wavefunction_1[i + 5] * wavefunction_1[i + 5];
        probabilities[i + 6] = wavefunction_1[i + 6] * wavefunction_1[i + 6];
        probabilities[i + 7] = wavefunction_1[i + 7] * wavefunction_1[i + 7]; 
        }
    for (i = size - (size % 8); i < (size); i++)
        probabilities[i] = wavefunction_1[i] * wavefunction_1[i];
    // end of vectorisation code    
    return(0);
    }
template <typename T>
T Slater_basis_set_calculations<T>::Probabilities_lenght(T* probabilities, unsigned int lenght_order,
unsigned int x_range, unsigned int y_range, unsigned int z_range, int x, int y, int z)
    { // Calculate average electron distance from point (x, y, z) for potential calculation
    unsigned int side, size;
    unsigned int i, j, k;
    unsigned int i_min, j_min, k_min;
    unsigned int i_max, j_max, k_max;
    T lenght;
    T l;
    T m;
    T n;
    
    i_min = lenght_order - z_range;
    j_min = lenght_order - y_range;
    k_min = lenght_order - x_range;
    i_max = lenght_order + z_range + 1;
    j_max = lenght_order + y_range + 1;
    k_max = lenght_order + x_range + 1;
    
    lenght = 0.00;
    side = (2 * lenght_order + 1);
    size = side * side * side;
    j = 0;
    k = 0;
    
    for (i = i_min; i < i_max; i++)
        {
        l = (i - lenght_order + z) * (i - lenght_order + z);
        for (j = j_min; j < j_max; j++)
            {
            m = (j - lenght_order + y) * (j - lenght_order + y);
            for (k = k_min; k < k_max; k++)
                {
                n = (k - lenght_order + x) * (k - lenght_order + x);
                if (l != 0 or m != 0 or n != 0)
                    {
                    lenght = lenght + (probabilities[(i * side * side) +
                    (j * side) + k]/sqrt(l + m + n));
                    }
                }
            }
        }
    lenght = 1/lenght * T(vector_lenght)/T(lenght_order) * Hartree_lenght; 
    return(lenght);
    }
template <typename T>
T Slater_basis_set_calculations<T>::Probabilities_lenght(T* probabilities, unsigned int lenght_order,
unsigned int x_range, unsigned int y_range, unsigned int z_range)
    { // Calculate electron distance from null point for potential calculation
    unsigned int side;
    unsigned int size;
    unsigned int i, j, k;
    unsigned int i_min, j_min, k_min;
    unsigned int i_max, j_max, k_max;
    unsigned int pre_index;
    T lenght;
    T* distances = results.lenghts[0];
    
    i_min = lenght_order - z_range;
    j_min = lenght_order - y_range;
    k_min = lenght_order - x_range;
    i_max = lenght_order + z_range + 1;
    j_max = lenght_order + y_range + 1;
    k_max = lenght_order + x_range + 1;
    
    lenght = 0.00;
    side = (2 * lenght_order + 1);
    size = side * side * side;
    j = 0;
    k = 0;
    
    for (i = i_min; i < i_max; i++)
        {
        for (j = j_min; j < j_max; j++)
            {
            for (k = k_min; k < k_max; k++)
                {
                pre_index = i * side * side + j * side;
                if (distances[pre_index + k] != 0)
                    {
                    lenght = lenght + (probabilities[(i * side * side) + (j * side) + k]/distances[pre_index + k]);
                    }
                }
            }
        }
    lenght = 1/lenght * T(vector_lenght)/T(lenght_order) * Hartree_lenght;
    return(lenght);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Probabilities_thread(T* Probabilities, unsigned int lenght_order,
unsigned int x_range, unsigned int y_range, unsigned int z_range, T* lenght)
    { // Calculate electron distance from null point for potential calculation for multithread calls
    lenght[0] = Probabilities_lenght(Probabilities, lenght_order, x_range, y_range, z_range);
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Gradient_thread(T* Gradient_1, T* wavefunction_2, unsigned int lenght_order)
    {
    unsigned int i, j;
    unsigned int side;
    unsigned int size;
    T x_gradient, y_gradient, z_gradient;
    T pixel_lenght;
    T Gradient;
    
    side = (2 * lenght_order + 1);
    size = side * side * side;
    pixel_lenght = T(vector_lenght)/T(lenght_order) * Hartree_lenght;
    
    for (i = 0; i < side * side; i++) // Optimalization for pipelining with non sequential reading of memory
        {
        if (i >= 1)
            x_gradient = (wavefunction_2[i + 1] - wavefunction_2[i - 1])/(pixel_lenght * 2);
        else
            x_gradient = (wavefunction_2[i + 1] - wavefunction_2[i])/pixel_lenght;
        
        if (i >= side)
            y_gradient = (wavefunction_2[i + side] - wavefunction_2[i - side])/(pixel_lenght * 2);
        else
            y_gradient = (wavefunction_2[i + side] - wavefunction_2[i])/pixel_lenght;
        
        z_gradient = (wavefunction_2[i + (side * side)] - wavefunction_2[i])/pixel_lenght;
        Gradient_1[i] = sqrt(((x_gradient * x_gradient) + (y_gradient * y_gradient) + (z_gradient * z_gradient))/3);
        }
    for (i = (side * side); i < size - (side * side); i++)
        {
        x_gradient = (wavefunction_2[i + 1] - wavefunction_2[i - 1])/(pixel_lenght * 2);
        y_gradient = (wavefunction_2[i + side] - wavefunction_2[i - side])/(pixel_lenght * 2);
        z_gradient = (wavefunction_2[i + (side * side)] - wavefunction_2[i - (side * side)])/(pixel_lenght * 2);
        Gradient_1[i] = sqrt(((x_gradient * x_gradient) + (y_gradient * y_gradient) + (z_gradient * z_gradient))/3);
        }
    for (i = size - (side * side); i < size ; i++)
        {
        if (size - i > 1)
            x_gradient = (wavefunction_2[i + 1] - wavefunction_2[i - 1])/(pixel_lenght * 2);
        else
            x_gradient = (wavefunction_2[i] - wavefunction_2[i - 1])/pixel_lenght;
        
        if (size - i > side)
            y_gradient = (wavefunction_2[i + side] - wavefunction_2[i - side])/(pixel_lenght * 2);
        else
            y_gradient = (wavefunction_2[i] - wavefunction_2[i - side])/pixel_lenght;
        
        z_gradient = (wavefunction_2[i] - wavefunction_2[i - (side * side)])/(pixel_lenght);
        Gradient_1[i] = sqrt(((x_gradient * x_gradient) + (y_gradient * y_gradient) + (z_gradient * z_gradient))/3);
        }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Integral_overlap(T* wavefunction_1, T* wavefunction_2, T* result,
unsigned int lenght_order, unsigned int x_range_1, unsigned int x_range_2, unsigned int y_range_1,
unsigned int y_range_2, unsigned int z_range_1, unsigned int z_range_2, T x, T y, T z) 
    { // Calculate the absolute value of overleap integral for multiplying the energy differences from fock matrix
    unsigned int i, j, k;
    unsigned int i_min, j_min, k_min;
    unsigned int i_max, j_max, k_max;
    unsigned int pre_index;
    unsigned int side;
    unsigned int size;
    
    unsigned int x_contraction, y_contraction, z_contraction;
    T* overlap_array;
    T overlap[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    T overlap_sum;
    
    x_contraction = abs(x * lenght_order/vector_lenght);
    y_contraction = abs(y * lenght_order/vector_lenght);
    z_contraction = abs(z * lenght_order/vector_lenght);
    side = (2 * lenght_order + 1);
    size = (side - x_contraction) * (side - y_contraction) * (side - z_contraction);
    if (x_contraction >= side or y_contraction >= side  or z_contraction >= side)
        return(-1);
    
    try {
        overlap_array = new T[size];
        }
    catch (int) {
        overlap_integral_exception_handle = 1;
        }
    if (x_contraction >= (lenght_order + 1) or y_contraction >= (lenght_order + 1) or z_contraction >= (lenght_order + 1))
        {
        result[0] = 0;
        return(-1);
        }
    // Detect indexes ranges for integration
    if (z_range_1 > z_range_2) {
        i_min = lenght_order - z_range_1;
        i_max = lenght_order + z_range_1 + 1;
        }
    else {
        i_min = lenght_order - z_range_2;
        i_max = lenght_order + z_range_2 + 1;
        }
    if (y_range_1 > y_range_2) {
        j_min = lenght_order - y_range_1;
        j_max = lenght_order + y_range_1 + 1;
        }
    else {
        j_min = lenght_order - y_range_2;
        j_max = lenght_order + y_range_2 + 1;
        }
    if (x_range_1 > x_range_2) {
        k_min = lenght_order - x_range_1;
        k_max = lenght_order + x_range_1 + 1;
        }
    else {
        k_min = lenght_order - x_range_2;
        k_max = lenght_order + x_range_2 + 1;
        }
    if (x == 0 and y == 0 and z == 0)
        Wavefunction_multiply(wavefunction_1, wavefunction_2, overlap_array, lenght_order);
    else
        Wavefunction_multiply(wavefunction_1, wavefunction_2, overlap_array, lenght_order, x, y, z);
    // vectorisation code
    if (x == 0 and y == 0 and z == 0)
        {
        for (i = i_min; i < i_max; i++)
            for (j = j_min; j < j_max; j++)
                {
                pre_index = i * side * side + j * side;
                for (k = k_min; k + 7 < k_max; k+=8)
                    {
                    overlap[0] += overlap_array[pre_index + k];
                    overlap[1] += overlap_array[pre_index + k + 1];
                    overlap[2] += overlap_array[pre_index + k + 2];
                    overlap[3] += overlap_array[pre_index + k + 3];
                    overlap[4] += overlap_array[pre_index + k + 4];
                    overlap[5] += overlap_array[pre_index + k + 5];
                    overlap[6] += overlap_array[pre_index + k + 6];
                    overlap[7] += overlap_array[pre_index + k + 7];
                    }
                for (k = k_max - (k_max % 8); k < k_max; k++)
                    overlap[0] += overlap_array[pre_index + k];
                }
        }
    else
        {
        for (i = 0; i + 7 < size; i = i + 8) // Multiplying of Wavefunction_1 and Wavefunction_1
            {
            overlap[0] += overlap_array[i];
            overlap[1] += overlap_array[i + 1];
            overlap[2] += overlap_array[i + 2];
            overlap[3] += overlap_array[i + 3];
            overlap[4] += overlap_array[i + 4];
            overlap[5] += overlap_array[i + 5];
            overlap[6] += overlap_array[i + 6];
            overlap[7] += overlap_array[i + 7]; 
            }
        for (i = size - (size % 8); i < (size); i++)
            overlap[0] = overlap[0] + overlap_array[i];
        }
    overlap_sum = overlap[0] + overlap[1] + overlap[2] + overlap[3] + overlap[4] + overlap[5] + overlap[6] + overlap[7];
    // end of vectorisation code    
    delete[] overlap_array;
    if ((not (isnan(overlap_sum))) and (not (isinf(overlap_sum)))) // Check for NaN and inf values
        result[0] = overlap_sum;
    
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Integrate_Integral_overlap(T* wavefunction_1, T* wavefunction_2, T* result,
unsigned int lenght_order, T x, T y, T z) 
    { // Calculate the absolute value of overleap integral for multiplying the energy differences from fock matrix
    unsigned int side;
    unsigned int size;
    unsigned int i, j, k, l, m, n;
    unsigned int x_contraction, y_contraction, z_contraction, x_side, y_side, z_side;
    int x_shift, y_shift, z_shift;
    T* overlap_array;
    T overlap[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    T overlap_sum;
    
    side = (2 * lenght_order + 1);
    x_shift = x * lenght_order/vector_lenght;
    y_shift = y * lenght_order/vector_lenght;
    z_shift = z * lenght_order/vector_lenght;
    x_contraction = abs(x_shift);
    y_contraction = abs(y_shift);
    z_contraction = abs(z_shift);
    x_side = side - x_contraction;
    y_side = side - y_contraction;
    z_side = side - z_contraction;
    
    T point_value_distance;
    T point_value;
    T average_lenght;
    T threshold = wavefunction_integration_threshold;
    
    size = (side - x_contraction) * (side - y_contraction) * (side - z_contraction);
    try {
        overlap_array = new T[size];
        }
    catch (int) {
        overlap_integral_exception_handle = 1;
        }
    if (x > (lenght_order + 1) or y > (lenght_order + 1) or z > (lenght_order + 1))
        {
        result[0] = 0;
        return(-1);
        }
        
    if (x == 0 and y == 0 and z == 0)
        Wavefunction_multiply(wavefunction_1, wavefunction_2, overlap_array, lenght_order);
    else
        Wavefunction_multiply(wavefunction_1, wavefunction_2, overlap_array, lenght_order, x, y, z);
    // vectorisation code    
    for (i = 0; i + 7 < size; i = i + 8) // Multiplying of Wavefunction_1 and Wavefunction_1
        {
        overlap[0] = overlap[0] + abs(overlap_array[i]);
        overlap[1] = overlap[1] + abs(overlap_array[i + 1]);
        overlap[2] = overlap[2] + abs(overlap_array[i + 2]);
        overlap[3] = overlap[3] + abs(overlap_array[i + 3]);
        overlap[4] = overlap[4] + abs(overlap_array[i + 4]);
        overlap[5] = overlap[5] + abs(overlap_array[i + 5]);
        overlap[6] = overlap[6] + abs(overlap_array[i + 6]);
        overlap[7] = overlap[7] + abs(overlap_array[i + 7]); 
        }
    for (i = size - (size % 8); i < (size); i++)
        overlap[0] = overlap[0] + abs(overlap_array[i]);
        
    overlap_sum = overlap[0] + overlap[1] + overlap[2] + overlap[3] + overlap[4] + overlap[5] + overlap[6] + overlap[7];
    // end of vectorisation code
    average_lenght = 0;
    for (k = 0; k < z_side; k++)
        for (j = 0; j < y_side; j++)
            for (i = 0; i < x_side; i++)
                if (abs(overlap_array[i + j * x_side + k * x_side * y_side]) > threshold)
                    for (n = 0; n < z_side; n++)
                        for (m = 0; m < y_side; m++)
                            for (l = 0; l < x_side; l++)
                                if (i != l or j != m or k != n)
                                    {
                                    point_value_distance = sqrt((int(i) - int(l)) * (int(i) - int(l)) +
                                    (int(j) - int(m)) * (int(j) - int(m)) + (int(k) - int(n)) * (int(k) - int(n)));
                                    point_value = 1/(point_value_distance * vector_lenght/lenght_order) *
                                    abs(overlap_array[i + j * x_side + k * x_side * y_side]) *
                                    abs(overlap_array[l + m * x_side + n * x_side * y_side]);
                                    average_lenght = average_lenght + point_value;
                                    }
    average_lenght = 1/average_lenght * (overlap_sum * overlap_sum) * Hartree_lenght;
    
    if ((not (isnan(average_lenght))) and (not (isinf(average_lenght)))) // Check for NaN and inf values
        result[0] = average_lenght;
    
    delete[] overlap_array;
    return(0);
    }
template <typename T>
inline int Slater_basis_set_calculations<T>::Integral_coulombic(T radius_1, T radius_2, T distance,  T* result, bool spin_bonded)
    { // Calculate the coulombic integral - potential energy of electron to proton of other atoms, aproximate
    T constant;
    T radius;
    T effective_lenght;
    
    constant = e*e/(4*Pi*E0);
    
    if (spin_bonded == true)
        {
        if (radius_1 >= radius_2)
            radius = radius_1 + (radius_2 * (Phi + 1)/3.00);
        else
            radius = radius_2 + (radius_1 * (Phi + 1)/3.00);
        }
    else
        {
        if (radius_1 >= radius_2)
            radius = radius_1 + (Phi - 1) * radius_2;
        else
            radius = radius_2 + (Phi - 1) * radius_1;
        }
    effective_lenght = sqrt((radius * radius) + (distance * distance));
    if ((not (isnan(constant/effective_lenght))) and (not (isinf(constant/effective_lenght)))) // Check for NaN and inf values
        result[0] = constant/effective_lenght;
    
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Integrate_Integral_coulombic(T* density_1, T* density_2, T* result,
unsigned int lenght_order, T x, T y, T z, unsigned int x_range_1, unsigned int x_range_2,
unsigned int y_range_1, unsigned int y_range_2, unsigned int z_range_1, unsigned int z_range_2) 
    { // Calculate the absolute value of overleap integral for multiplying the energy differences from fock matrix
    unsigned int side;
    unsigned int i, j, k, l, m, n;
    unsigned int i_min, j_min, k_min, l_min, m_min, n_min;
    unsigned int i_max, j_max, k_max, l_max, m_max, n_max;
    unsigned int pre_index, pre_index_2;
    int x_shift, y_shift, z_shift;
    int y_2_z_2;
    
    T constant;
    T radius;
    T point_value_distance;
    T point_value;
    T average_lenght;
    T average_lenghts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    T* reverse_relative_lenghts = nullptr;
    
    if (small_results.relative_lenghts.size() > 0 and small_results.lenght_orders.size() > 0)
        if (small_results.lenght_orders[0] == lenght_order)
            reverse_relative_lenghts = small_results.relative_lenghts[0];
    
    side = (2 * lenght_order + 1);
    x_shift = x * lenght_order/vector_lenght;
    y_shift = y * lenght_order/vector_lenght;
    z_shift = z * lenght_order/vector_lenght;
    
    
    i_max = lenght_order + x_range_1 + 1;
    j_max = lenght_order + y_range_1 + 1;
    k_max = lenght_order + z_range_1 + 1;
    l_max = lenght_order + x_range_2 + 1;
    m_max = lenght_order + y_range_2 + 1;
    n_max = lenght_order + z_range_2 + 1;
    i_min = lenght_order - x_range_1;
    j_min = lenght_order - y_range_1;
    k_min = lenght_order - z_range_1;
    l_min = lenght_order - x_range_2;
    m_min = lenght_order - y_range_2;
    n_min = lenght_order - z_range_2;
    
    // integration
    // vectorisation code
    if (x == 0 and y == 0 and z == 0 and reverse_relative_lenghts != nullptr)
        {
        // optimized integration with zero of coordinate distance
        for (k = k_min; k < k_max; k++)
            for (j = j_min; j < j_max; j++)
                for (i = i_min; i < i_max; i++)
                    for (n = n_min; n < n_max; n++)
                        {
                        pre_index = k * side * side * side * side * side + j * side * side * side * side +
                        i * side * side * side + n * side * side;
                        for (m = m_min; m < m_max; m++)
                            {
                            pre_index_2 = m * side;
                            for (l = l_min; l + 7 < l_max; l+=8)
                                {
                                average_lenghts[0] += reverse_relative_lenghts[pre_index + pre_index_2 + l] *
                                density_1[i + j * side + k * side * side] *
                                density_2[l + m * side + n * side * side];
                                average_lenghts[1] += reverse_relative_lenghts[pre_index + pre_index_2 + l + 1] *
                                density_1[i + j * side + k * side * side] *
                                density_2[l + m * side + n * side * side + 1];
                                average_lenghts[2] += reverse_relative_lenghts[pre_index + pre_index_2 + l + 2] *
                                density_1[i + j * side + k * side * side] *
                                density_2[l + m * side + n * side * side + 2];
                                average_lenghts[3] += reverse_relative_lenghts[pre_index + pre_index_2 + l + 3] *
                                density_1[i + j * side + k * side * side] *
                                density_2[l + m * side + n * side * side + 3];
                                average_lenghts[4] += reverse_relative_lenghts[pre_index + pre_index_2 + l + 4] *
                                density_1[i + j * side + k * side * side] *
                                density_2[l + m * side + n * side * side + 4];
                                average_lenghts[5] += reverse_relative_lenghts[pre_index + pre_index_2 + l + 5] *
                                density_1[i + j * side + k * side * side] *
                                density_2[l + m * side + n * side * side + 5];
                                average_lenghts[6] += reverse_relative_lenghts[pre_index + pre_index_2 + l + 6] *
                                density_1[i + j * side + k * side * side] *
                                density_2[l + m * side + n * side * side + 6];
                                average_lenghts[7] += reverse_relative_lenghts[pre_index + pre_index_2 + l + 7] *
                                density_1[i + j * side + k * side * side] *
                                density_2[l + m * side + n * side * side + 7];
                                }
                            for (l = l_max - (l_max % 8); l < l_max; l++)
                                {
                                average_lenghts[0] += reverse_relative_lenghts[pre_index + pre_index_2 + l] *
                                density_1[i + j * side + k * side * side] *
                                density_2[l + m * side + n * side * side];
                                }
                            }
                        }
        average_lenght = 1/(average_lenghts[0] + average_lenghts[1] + average_lenghts[2] + average_lenghts[3] +
        average_lenghts[4] + average_lenghts[5] + average_lenghts[6] + average_lenghts[7])
        *  T(lenght_order)/T(vector_lenght)  * Hartree_lenght;
        }
    else
        { // standard integration with non-zero of coordinate distance
        for (k = k_min; k < k_max; k++)
            for (j = j_min; j < j_max; j++)
                for (i = i_min; i < i_max; i++)
                    for (n = n_min; n < n_max; n++)
                        {
                        y_2_z_2 = (int(j) - int(m) - int(y_shift)) * (int(j) - int(m) - int(y_shift)) +
                        (int(k) - int(n) - int(z_shift)) * (int(k) - int(n) - int(z_shift));
                        for (m = m_min; m < m_max; m++)
                            {
                            for (l = l_min; l < l_max; l++)
                                {
                                point_value_distance = sqrt((int(i) - int(l) - int(x_shift)) *
                                (int(i) - int(l) - int(x_shift)) + y_2_z_2);
                                point_value = 1/(point_value_distance) *
                                density_1[i + j * side + k * side * side] *
                                density_2[l + m * side + n * side * side];
                                average_lenghts[0] += point_value;
                                }
                            }
                        }
        average_lenght = 1/(average_lenghts[0]) *  T(lenght_order)/T(vector_lenght)  * Hartree_lenght;
        }
    // End of vectorisation code
    constant = e*e/(4*Pi*E0);
    radius = average_lenght;
    
    if ((not (isnan(constant/radius))) and (not (isinf(constant/radius)))) // Check for NaN and inf values
        result[0] = constant/radius;
    
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Integral_nucleus_atraction(T probabilities_lenght, T multiplier,  T* result, T* lenght, unsigned int Z)
    { // Calculate the nucleus attraction integral - potential energy of electron to atom nucleus
    int lenght_X, lenght_Y, lenght_Z; // Z is equal the number of protons - number of eletrons without coulombic integral
    T constant;
    T radius;
    
    constant = -e * e/(4 * Pi * E0);
    radius = probabilities_lenght/multiplier;
    lenght[0] = radius;
    
    if ((not (isnan(constant * Z/radius))) and (not (isinf(constant * Z/radius))))
        result[0] = result[0] + constant * Z/radius; // Check for NaN and inf values
    
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Integrate_Integral_nucleus_atraction(T* probabilities,
T* result, T* lenght, unsigned int lenght_order, unsigned int x_range, unsigned int y_range, unsigned int z_range,
T lenght_x, T lenght_y, T lenght_z, unsigned int Z)
    { // Calculate the nucleus attraction integral - potential energy of electron to atom nucleus
    int lenght_X, lenght_Y, lenght_Z; // Z is equal the number of protons - number of eletrons without coulombic integral
    T constant;
    T radius;
    
    constant = -e * e/(4 * Pi * E0);
    lenght_X = lenght_x * lenght_order/vector_lenght;
    lenght_Y = lenght_y * lenght_order/vector_lenght;
    lenght_Z = lenght_z * lenght_order/vector_lenght; 
    radius = Probabilities_lenght(probabilities, lenght_order, x_range, y_range, z_range, lenght_X, lenght_Y, lenght_Z);
    lenght[0] = radius;
        
    if ((not (isnan(constant * Z/radius))) and (not (isinf(constant * Z/radius))))
        result[0] = result[0] + constant * Z/radius; // Check for NaN and inf values
    
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Integral_kinetic(T* Gradient_1, T* Gradient_2, T* result,
unsigned int lenght_order, T d_x, T d_y, T d_z, unsigned int x_range_1, unsigned int x_range_2,
unsigned int y_range_1, unsigned int y_range_2, unsigned int z_range_1, unsigned int z_range_2)
    {
    unsigned int i, j, k;
    unsigned int i_min, j_min, k_min;
    unsigned int i_max, j_max, k_max;
    unsigned int x_threshold_1, x_threshold_2;
    unsigned int y_threshold_1, y_threshold_2;
    unsigned int z_threshold_1, z_threshold_2;
    
    unsigned int side;
    unsigned int x_side, y_side, z_side;
    unsigned int x_contraction, y_contraction, z_contraction;
    unsigned int x_1_min, y_1_min, z_1_min;
    unsigned int x_2_min, y_2_min, z_2_min;
    int x, y, z;
    
    T result_array[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    
    x = d_x * lenght_order/vector_lenght;
    y = d_y * lenght_order/vector_lenght;
    z = d_z * lenght_order/vector_lenght;
    side = (2 * lenght_order + 1);
    x_contraction = abs(x);
    y_contraction = abs(y);
    z_contraction = abs(z);
    if (x_contraction > side)
        x_contraction = side;
    if (y_contraction > side)
        y_contraction = side;
    if (z_contraction > side)
        z_contraction = side;
    
    x_side = side - x_contraction;
    y_side = side - y_contraction;
    z_side = side - z_contraction;
    
    if (x == 0 and y == 0 and z == 0)
        {
        if (x_range_1 > x_range_2) {
            i_min = lenght_order - z_range_1;
            i_max = lenght_order + z_range_1 + 1;}
        else {
            i_min = lenght_order - z_range_2;
            i_max = lenght_order + z_range_2 + 1;}
    
        if (y_range_1 > y_range_2) {
           j_min = lenght_order - y_range_1;
            j_max = lenght_order + y_range_1 + 1;}
        else {
            j_min = lenght_order - y_range_2;
            j_max = lenght_order + y_range_2 + 1;}
    
        if (z_range_1 > z_range_2) {
            k_min = lenght_order - x_range_1;
            k_max = lenght_order + x_range_1 + 1;}
        else {
            k_min = lenght_order - x_range_2;
            k_max = lenght_order + x_range_2 + 1;}
        }
    else
        {
        if (x + x_range_1 > x_range_2) {
            k_min = lenght_order - x_range_2;}
        else {
            k_min = lenght_order - x_range_1;}
        
        if (y + y_range_1 > y_range_2) {
            j_min = lenght_order - y_range_2;}
        else {
            j_min = lenght_order - y_range_1;}
        
        if (z + z_range_1 > z_range_2) {
            i_min = lenght_order - z_range_2;}
        else {
            i_min = lenght_order - z_range_1;}
            
        if (x + x_range_2 > x_range_1) {
            k_max = lenght_order + x_range_1;}
        else {
            k_max = lenght_order + x_range_2;}
        
        if (y + y_range_2 > y_range_1) {
            j_max = lenght_order + y_range_1;}
        else {
            j_max = lenght_order + y_range_2;}
        
        if (z + z_range_2 > z_range_1) {
            i_max = lenght_order + z_range_1;}
        else {
            i_max = lenght_order + z_range_2;}
        }
    if (x_contraction >= side or y_contraction >= side or z_contraction >= side)
        return(-1);
    
    
    if (x < 0)
        {
        x_1_min = 0;
        x_2_min = x_contraction;
        }
    else
        {
        x_1_min = x_contraction;
        x_2_min = 0;
        }
    if (y < 0)
        {
        y_1_min = 0;
        y_2_min = y_contraction;
        }
    else
        {
        y_1_min = y_contraction;
        y_2_min = 0;
        }
    if (z < 0)
        {
        z_1_min = 0;
        z_2_min = z_contraction;
        }
    else
        {
        z_1_min = z_contraction;
        z_2_min = 0;
        }
    // vectorisation code
    if (d_x != 0 or d_y != 0 or d_z != 0)
        for (i = i_min; i + z_contraction < (1 + i_max); i++)
            {
            for (j = j_min; j + y_contraction < (1 + j_max); j++)
                {
                for (k = k_min; k + x_contraction < (1 + k_max) + 7; k += 8)
                    {
                    result_array[0] = result_array[0] +
                    Gradient_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min)] *
                    Gradient_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min)];
                    result_array[1] = result_array[1] +
                    Gradient_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 1] *
                    Gradient_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 1];
                    result_array[2] = result_array[2] + 
                    Gradient_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 2] *
                    Gradient_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 2];
                    result_array[3] = result_array[3] + 
                    Gradient_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 3] *
                    Gradient_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 3];
                    result_array[4] = result_array[4] + 
                    Gradient_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 4] *
                    Gradient_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 4];
                    result_array[5] = result_array[5] + 
                    Gradient_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 5] *
                    Gradient_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 5];
                    result_array[6] = result_array[6] + 
                    Gradient_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 6] *
                    Gradient_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 6];
                    result_array[7] = result_array[7] + 
                    Gradient_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min) + 7] *
                    Gradient_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min) + 7];
                    }
                if (x_contraction < (1 + j_max))
                    for (k = 1 + k_max - x_contraction - ((1 + k_max - x_contraction) % 8);
                    k < 1 + k_max - x_contraction; k++)
                        {
                        result_array[0] = result_array[0] + 
                        Gradient_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min)] *
                        Gradient_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min)];
                        }
                }
            }
    else
        for (i = i_min; i < i_max; i++)
            {
            for (j = j_min; j < j_max; j++)
                {
                for (k = k_min; k + 7 < k_max; k += 8)
                    {
                    result_array[0] = result_array[0] + Gradient_1[(i * side * side) + (j * side) + k] *
                    Gradient_2[(i * side * side) + (j * side) + k];
                    result_array[1] = result_array[1] + Gradient_1[(i * side * side) + (j * side) + k + 1] *
                    Gradient_2[(i * side * side) + (j * side) + k + 1];
                    result_array[2] = result_array[2] + Gradient_1[(i * side * side) + (j * side) + k + 2] *
                    Gradient_2[(i * side * side) + (j * side) + k + 2];
                    result_array[3] = result_array[3] + Gradient_1[(i * side * side) + (j * side) + k + 3] *
                    Gradient_2[(i * side * side) + (j * side) + k + 3];
                    result_array[4] = result_array[3] + Gradient_1[(i * side * side) + (j * side) + k + 4] *
                    Gradient_2[(i * side * side) + (j * side) + k + 4];
                    result_array[5] = result_array[5] + Gradient_1[(i * side * side) + (j * side) + k + 5] *
                    Gradient_2[(i * side * side) + (j * side) + k + 5];
                    result_array[6] = result_array[6] + Gradient_1[(i * side * side) + (j * side) + k + 6] *
                    Gradient_2[(i * side * side) + (j * side) + k + 6];
                    result_array[7] = result_array[7] + Gradient_1[(i * side * side) + (j * side) + k + 7] *
                    Gradient_2[(i * side * side) + (j * side) + k + 7];
                    }
                for (k = k_max - (k_max % 8); k < k_max; k++)
                    {
                    result_array[0] = result_array[0] + Gradient_1[(i * side * side) + (j * side) + k] *
                    Gradient_2[(i * side * side) + (j * side) + k];
                    }
                }
            }
    result[0] = (result_array[0] + result_array[1] + result_array[2] + result_array[3] + result_array[4] +
    result_array[5] + result_array[6] + result_array[7]) * (-h * h)/(8 * Pi * Pi * me);
        // end of vectorisation code
    return(0);
    }
template <typename T>
T Slater_basis_set_calculations<T>::Rydberg_energy(unsigned int Z, unsigned int n)
    {
    T energy;
    
    energy = -(Z * Z * Hartree_energy_constant/(2 * n * n));
    return(energy);
    }
template <typename T>
T Slater_basis_set_calculations<T>::Spin_moment_energy(T s, T B0)
    {
    T energy;
    T S;
    
    if (s * (s + 1) != 0)
        S = (s * (s + 1)) * h/(4 * Pi);
    else
        S = 0;
        
    energy = (S * e * h)/(4 * Pi * me) * 2.00232 * B0; // including Bohr magnetron
    return(energy);
    }
template <typename T>
T Slater_basis_set_calculations<T>::Orbital_magnetic_field(T potential_energy, T radius, int l)
    {
    T B;
    T L;
    
    L = l * (l + 1);
    B = potential_energy/(radius * radius * me * e * c * c) * L;
    return(B);
    }
// End of Section 3 - mathematical operations for wavefunctions, probabiliti densities and integrals
//Section 4 - generating lists of electrons
template <typename T>
int Slater_basis_set_calculations<T>::Quantum_numbers_to_orbitals(unsigned int n, unsigned int l, int fulness,
atom_orbitals* atom_orbitals_PTR)
    { // according to Hund rule, there are not control of fulness parameter overload
    if (l == 0) // s orbitals
        {
        atom_orbitals_PTR->n.push_back(n);
        atom_orbitals_PTR->l.push_back(l);
        atom_orbitals_PTR->m.push_back(0);
        atom_orbitals_PTR->s.push_back(0.5);
        atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
        atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
        atom_orbitals_PTR->bonding.push_back(-1);
        atom_orbitals_PTR->paired_with_previous.push_back(false);
        if (fulness == 2)
            {
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(0);
            atom_orbitals_PTR->s.push_back(-0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(true);
            }
        }
    if (l == 1) // p orbitals
        {
        if (fulness >= 1)
            { // px orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(1);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/6) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(1);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        if (fulness >= 2)
            { // pz orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(0);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/5) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(0);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            
            }
        if (fulness >= 3)
            { // py orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(-1);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false); 
            if ((fulness/4) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(-1);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        }
    if (l == 2) // d orbitals
        {
        if (fulness >= 1)
            { // dx2-y2 orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(2);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/10) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(2);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        if (fulness >= 2)
            { // dxz orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(1);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/9) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(1);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        if (fulness >= 3)
            { // dz2 orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(0);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/8) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(0);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        if (fulness >= 4)
            {// dyz orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(-1);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false); 
            if ((fulness/7) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(-1);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        if (fulness >= 5)
            { // dxy orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(-2);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false); 
            if ((fulness/6) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(-2);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        }
    if (l == 3) // f orbitals
        {
        if (fulness >= 1)
            { // fx(x2-3y2) orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(3);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/14) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(3);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        if (fulness >= 2)
            {  // fz(x2-y2) orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(2);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/13) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(2);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        if (fulness >= 3)
            { // fxz2 orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(1);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/12) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(1);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        if (fulness >= 4)
            { // fz3 orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(0);
            atom_orbitals_PTR->s.push_back(0.5); 
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/11) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(0);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        if (fulness >= 5)
            { // fyz2 orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(-1);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/10) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(-1);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        if (fulness >= 6)
            { // fz(xy) orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(-2);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/9) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(-2);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        if (fulness >= 7)
            { // fy(x2-3y2) orbital
            atom_orbitals_PTR->n.push_back(n);
            atom_orbitals_PTR->l.push_back(l);
            atom_orbitals_PTR->m.push_back(-3);
            atom_orbitals_PTR->s.push_back(0.5);
            atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
            atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
            atom_orbitals_PTR->bonding.push_back(-1);
            atom_orbitals_PTR->paired_with_previous.push_back(false);
            if ((fulness/8) == 1)
                {
                atom_orbitals_PTR->n.push_back(n);
                atom_orbitals_PTR->l.push_back(l);
                atom_orbitals_PTR->m.push_back(-3);
                atom_orbitals_PTR->s.push_back(-0.5);
                atom_orbitals_PTR->wavefunction_lenght_multipliers.push_back(1);
                atom_orbitals_PTR->wavefunction_coefficients.push_back(1);
                atom_orbitals_PTR->bonding.push_back(-1);
                atom_orbitals_PTR->paired_with_previous.push_back(true);
                }
            }
        }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Atoms_to_valence_orbitals(string atom, atom_orbitals* atom_orbitals_PTR)
    {
    if (atom == "H") // 1s1
        {
        Quantum_numbers_to_orbitals(1, 0, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 1;
        atom_orbitals_PTR->reduced_Z = 1;
        atom_orbitals_PTR->electronegativity = 2.2;
        }
    if (atom == "He") // 1s2
        {
        Quantum_numbers_to_orbitals(1, 0, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 2;
        atom_orbitals_PTR->reduced_Z = 2;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Li") // 2s1
        {
        Quantum_numbers_to_orbitals(2, 0, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 3;
        atom_orbitals_PTR->reduced_Z = 1;
        atom_orbitals_PTR->electronegativity = 0.98;
        }
    if (atom == "Be") // 2s2
        {
        Quantum_numbers_to_orbitals(2, 0, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 4;
        atom_orbitals_PTR->reduced_Z = 2;
        atom_orbitals_PTR->electronegativity = 1.57;
        }
    if (atom == "B") // 2s2 2p1
        {
        Quantum_numbers_to_orbitals(2, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(2, 1, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 5;
        atom_orbitals_PTR->reduced_Z = 3;
        atom_orbitals_PTR->electronegativity = 2.04;
        }
    if (atom == "C") // 2s2 2p2
        {
        Quantum_numbers_to_orbitals(2, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(2, 1, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 6;
        atom_orbitals_PTR->reduced_Z = 4;
        atom_orbitals_PTR->electronegativity = 2.55;
        }
    if (atom == "N") // 2s2 2p3
        {
        Quantum_numbers_to_orbitals(2, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(2, 1, 3, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 7;
        atom_orbitals_PTR->reduced_Z = 5;
        atom_orbitals_PTR->electronegativity = 3.04;
        }
    if (atom == "O") // 2s2 2p4
        {
        Quantum_numbers_to_orbitals(2, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(2, 1, 4, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 8;
        atom_orbitals_PTR->reduced_Z = 6;
        atom_orbitals_PTR->electronegativity = 3.44;
        }
    if (atom == "F") // 2s2 2p5
        {
        Quantum_numbers_to_orbitals(2, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(2, 1, 5, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 9;
        atom_orbitals_PTR->reduced_Z = 7;
        atom_orbitals_PTR->electronegativity = 3.98;
        }
    if (atom == "Ne") // 2s2 2p6
        {
        Quantum_numbers_to_orbitals(2, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(2, 1, 6, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 10;
        atom_orbitals_PTR->reduced_Z = 8;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Na") // 3s1
        {
        Quantum_numbers_to_orbitals(3, 0, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 11;
        atom_orbitals_PTR->reduced_Z = 1;
        atom_orbitals_PTR->electronegativity = 0.93;
        }
    if (atom == "Mg") // 3s2
        {
        Quantum_numbers_to_orbitals(3, 0, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 12;
        atom_orbitals_PTR->reduced_Z = 2;
        atom_orbitals_PTR->electronegativity = 1.31;
        }
    if (atom == "Al") // 3s2 3p1
        {
        Quantum_numbers_to_orbitals(3, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 1, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 13;
        atom_orbitals_PTR->reduced_Z = 3;
        atom_orbitals_PTR->electronegativity = 1.61;
        }
    if (atom == "Si") // 3s2 3p2
        {
        Quantum_numbers_to_orbitals(3, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 1, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 14;
        atom_orbitals_PTR->reduced_Z = 4;
        atom_orbitals_PTR->electronegativity = 1.9;
        }
    if (atom == "P") // 3s2 3p3
        {
        Quantum_numbers_to_orbitals(3, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 1, 3, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 15;
        atom_orbitals_PTR->reduced_Z = 5;
        atom_orbitals_PTR->electronegativity = 2.19;
        }
    if (atom == "S") // 3s2 3p4
        {
        Quantum_numbers_to_orbitals(3, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 1, 4, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 16;
        atom_orbitals_PTR->reduced_Z = 6;
        atom_orbitals_PTR->electronegativity = 2.58;
        }
    if (atom == "Cl") // 3s2 3p5
        {
        Quantum_numbers_to_orbitals(3, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 1, 5, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 17;
        atom_orbitals_PTR->reduced_Z = 7;
        atom_orbitals_PTR->electronegativity = 3.16;
        }
    if (atom == "Ar") // 3s2 3p6
        {
        Quantum_numbers_to_orbitals(3, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 1, 6, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 18;
        atom_orbitals_PTR->reduced_Z = 8;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "K") // 4s1
        {
        Quantum_numbers_to_orbitals(4, 0, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 19;
        atom_orbitals_PTR->reduced_Z = 1;
        atom_orbitals_PTR->electronegativity = 0.82;
        }
    if (atom == "Ca") // 4s2
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 20;
        atom_orbitals_PTR->reduced_Z = 2;
        atom_orbitals_PTR->electronegativity = 1;
        }
    if (atom == "Sc") // 4s2 3d1
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 21;
        atom_orbitals_PTR->reduced_Z = 3;
        atom_orbitals_PTR->electronegativity = 1.36;
        }
    if (atom == "Ti") // 4s2 3d2
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 22;
        atom_orbitals_PTR->reduced_Z = 4;
        atom_orbitals_PTR->electronegativity = 1.54;
        }
    if (atom == "V") // 4s2 3d3
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 3, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 23;
        atom_orbitals_PTR->reduced_Z = 5;
        atom_orbitals_PTR->electronegativity = 1.63;
        }
    if (atom == "Cr") // 4s1 3d5
        {
        Quantum_numbers_to_orbitals(4, 0, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 5, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 24;
        atom_orbitals_PTR->reduced_Z = 6;
        atom_orbitals_PTR->electronegativity = 1.66;
        }
    if (atom == "Mn") // 4s2 3d5
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 5, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 25;
        atom_orbitals_PTR->reduced_Z = 7;
        atom_orbitals_PTR->electronegativity = 1.55;
        }
    if (atom == "Fe") // 4s2 3d6
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 6, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 26;
        atom_orbitals_PTR->reduced_Z = 8;
        atom_orbitals_PTR->electronegativity = 1.83;
        }
    if (atom == "Co") // 4s2 3d7
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 7, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 27;
        atom_orbitals_PTR->reduced_Z = 9;
        atom_orbitals_PTR->electronegativity = 1.88;
        }
    if (atom == "Ni") // 4s2 3d8
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 8, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 28;
        atom_orbitals_PTR->reduced_Z = 10;
        atom_orbitals_PTR->electronegativity = 1.91;
        }
    if (atom == "Cu") // 4s1 3d10
        {
        Quantum_numbers_to_orbitals(4, 0, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 29;
        atom_orbitals_PTR->reduced_Z = 11;
        atom_orbitals_PTR->electronegativity = 1.9;
        }
    if (atom == "Zn") // 4s2 3d10
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 30;
        atom_orbitals_PTR->reduced_Z = 12;
        atom_orbitals_PTR->electronegativity = 1.65;
        }
    if (atom == "Ga") // 4s2 3d10 4p1
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 1, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 31;
        atom_orbitals_PTR->reduced_Z = 13;
        atom_orbitals_PTR->electronegativity = 1.81;
        }
    if (atom == "Ge") // 4s2 3d10 4p2
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 1, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 32;
        atom_orbitals_PTR->reduced_Z = 14;
        atom_orbitals_PTR->electronegativity = 2.01;
        }
    if (atom == "As") // 4s2 3d10 4p3
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 1, 3, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 33;
        atom_orbitals_PTR->reduced_Z = 15;
        atom_orbitals_PTR->electronegativity = 2.18;
        }
    if (atom == "Se") // 4s2 3d10 4p4
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 1, 4, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 34;
        atom_orbitals_PTR->reduced_Z = 16;
        atom_orbitals_PTR->electronegativity = 2.55;
        }
    if (atom == "Br") // 4s2 3d10 4p5
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 1, 5, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 35;
        atom_orbitals_PTR->reduced_Z = 17;
        atom_orbitals_PTR->electronegativity = 2.96;
        }
    if (atom == "Kr") // 4s2 3d10 4p6
        {
        Quantum_numbers_to_orbitals(4, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(3, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 1, 6, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 36;
        atom_orbitals_PTR->reduced_Z = 18;
        atom_orbitals_PTR->electronegativity = 3;
        }
    if (atom == "Rb") // 5s1
        {
        Quantum_numbers_to_orbitals(5, 0, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 37;
        atom_orbitals_PTR->reduced_Z = 1;
        atom_orbitals_PTR->electronegativity = 0.82;
        }
    if (atom == "Sr") // 5s2
        {
        Quantum_numbers_to_orbitals(5, 0, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 38;
        atom_orbitals_PTR->reduced_Z = 2;
        atom_orbitals_PTR->electronegativity = 0.95;
        }
    if (atom == "Y") // 5s2 4d1
        {
        Quantum_numbers_to_orbitals(5, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 39;
        atom_orbitals_PTR->reduced_Z = 3;
        atom_orbitals_PTR->electronegativity = 1.22;
        }
    if (atom == "Zr") // 5s2 4d2
        {
        Quantum_numbers_to_orbitals(5, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 40;
        atom_orbitals_PTR->reduced_Z = 4;
        atom_orbitals_PTR->electronegativity = 1.33;
        }
    if (atom == "Nb") // 5s2 4d3
        {
        Quantum_numbers_to_orbitals(5, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 3, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 41;
        atom_orbitals_PTR->reduced_Z = 5;
        atom_orbitals_PTR->electronegativity = 1.6;
        }
    if (atom == "Mo") // 5s1 4d5
        {
        Quantum_numbers_to_orbitals(5, 0, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 5, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 42;
        atom_orbitals_PTR->reduced_Z = 6;
        atom_orbitals_PTR->electronegativity = 2.16;
        }
    if (atom == "Tc") // 5s1 4d6
        {
        Quantum_numbers_to_orbitals(5, 0, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 6, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 43;
        atom_orbitals_PTR->reduced_Z = 7;
        atom_orbitals_PTR->electronegativity = 1.9;
        }
    if (atom == "Ru") // 5s1 4d7
        {
        Quantum_numbers_to_orbitals(5, 0, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 7, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 44;
        atom_orbitals_PTR->reduced_Z = 8;
        atom_orbitals_PTR->electronegativity = 2.2;
        }
    if (atom == "Rh") // 5s1 4d8
        {
        Quantum_numbers_to_orbitals(5, 0, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 8, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 45;
        atom_orbitals_PTR->reduced_Z = 9;
        atom_orbitals_PTR->electronegativity = 2.28;
        }
    if (atom == "Pd") // 4d10
        {
        Quantum_numbers_to_orbitals(4, 2, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 46;
        atom_orbitals_PTR->reduced_Z = 10;
        atom_orbitals_PTR->electronegativity = 2.2;
        }
    if (atom == "Ag") // 5s1 4d10
        {
        Quantum_numbers_to_orbitals(5, 0, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 47;
        atom_orbitals_PTR->reduced_Z = 11;
        atom_orbitals_PTR->electronegativity = 1.93;
        }
    if (atom == "Cd") // 5s1 4d10
        {
        Quantum_numbers_to_orbitals(5, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 47;
        atom_orbitals_PTR->reduced_Z = 11;
        atom_orbitals_PTR->electronegativity = 1.93;
        }
    if (atom == "In") // 5s2 5p1 4d10
        {
        Quantum_numbers_to_orbitals(5, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 1, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 49;
        atom_orbitals_PTR->reduced_Z = 13;
        atom_orbitals_PTR->electronegativity = 1.78;
        }
    if (atom == "Sn") // 5s2 5p2 4d10
        {
        Quantum_numbers_to_orbitals(5, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 1, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 50;
        atom_orbitals_PTR->reduced_Z = 14;
        atom_orbitals_PTR->electronegativity = 1.96;
        }
    if (atom == "Sb") // 5s2 5p3 4d10
        {
        Quantum_numbers_to_orbitals(5, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 1, 3, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 51;
        atom_orbitals_PTR->reduced_Z = 15;
        atom_orbitals_PTR->electronegativity = 2.05;
        }
    if (atom == "Te") // 5s2 5p4 4d10
        {
        Quantum_numbers_to_orbitals(5, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 1, 4, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 52;
        atom_orbitals_PTR->reduced_Z = 16;
        atom_orbitals_PTR->electronegativity = 2.1;
        }
    if (atom == "I") // 5s2 5p5 4d10
        {
        Quantum_numbers_to_orbitals(5, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 1, 5, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 53;
        atom_orbitals_PTR->reduced_Z = 17;
        atom_orbitals_PTR->electronegativity = 2.66;
        }
    if (atom == "Xe") // 5s2 5p6 4d10
        {
        Quantum_numbers_to_orbitals(5, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 1, 6, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 2, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 54;
        atom_orbitals_PTR->reduced_Z = 18;
        atom_orbitals_PTR->electronegativity = 2.6;
        }
    if (atom == "Cs") // 6s1
        {
        Quantum_numbers_to_orbitals(6, 0, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 55;
        atom_orbitals_PTR->reduced_Z = 1;
        atom_orbitals_PTR->electronegativity = 0.79;
        }
    if (atom == "Ba") // 6s2
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 56;
        atom_orbitals_PTR->reduced_Z = 2;
        atom_orbitals_PTR->electronegativity = 0.89;
        }
    if (atom == "La") // 6s2 5d1
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 57;
        atom_orbitals_PTR->reduced_Z = 3;
        atom_orbitals_PTR->electronegativity = 1.1;
        }
    if (atom == "Ce") // 6s2 5d1 4f1
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 58;
        atom_orbitals_PTR->reduced_Z = 4;
        atom_orbitals_PTR->electronegativity = 1.12;
        }
    if (atom == "Pr") // 6s2 4f3
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 3, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 59;
        atom_orbitals_PTR->reduced_Z = 5;
        atom_orbitals_PTR->electronegativity = 1.13;
        }
    if (atom == "Nd") // 6s2 4f4
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 4, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 60;
        atom_orbitals_PTR->reduced_Z = 6;
        atom_orbitals_PTR->electronegativity = 1.14;
        }
    if (atom == "Pm") // 6s2 4f5
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 5, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 61;
        atom_orbitals_PTR->reduced_Z = 7;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Sm") // 6s2 4f6
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 6, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 62;
        atom_orbitals_PTR->reduced_Z = 8;
        atom_orbitals_PTR->electronegativity = 1.17;
        }
    if (atom == "Eu") // 6s2 4f7
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 7, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 63;
        atom_orbitals_PTR->reduced_Z = 9;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Gd") // 6s2 5d1 4f7
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 7, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 64;
        atom_orbitals_PTR->reduced_Z = 10;
        atom_orbitals_PTR->electronegativity = 1.2;
        }
    if (atom == "Tb") // 6s2 4f9
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 9, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 65;
        atom_orbitals_PTR->reduced_Z = 11;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Dy") // 6s2 4f10
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 66;
        atom_orbitals_PTR->reduced_Z = 12;
        atom_orbitals_PTR->electronegativity = 1.22;
        }
    if (atom == "Ho") // 6s2 4f11
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 11, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 67;
        atom_orbitals_PTR->reduced_Z = 13;
        atom_orbitals_PTR->electronegativity = 1.23;
        }
    if (atom == "Er") // 6s2 4f12
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 12, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 68;
        atom_orbitals_PTR->reduced_Z = 14;
        atom_orbitals_PTR->electronegativity = 1.24;
        }
    if (atom == "Tm") // 6s2 4f13
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 13, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 69;
        atom_orbitals_PTR->reduced_Z = 15;
        atom_orbitals_PTR->electronegativity = 1.25;
        }
    if (atom == "Yb") // 6s2 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 70;
        atom_orbitals_PTR->reduced_Z = 16;
        atom_orbitals_PTR->electronegativity = 1.1;
        }
    if (atom == "Lu") // 6s2 5d1 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 71;
        atom_orbitals_PTR->reduced_Z = 17;
        atom_orbitals_PTR->electronegativity = 1.27;
        }
    if (atom == "Hf") // 6s2 5d2 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 72;
        atom_orbitals_PTR->reduced_Z = 18;
        atom_orbitals_PTR->electronegativity = 1.3;
        }
    if (atom == "Ta") // 6s2 5d3 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 3, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 73;
        atom_orbitals_PTR->reduced_Z = 19;
        atom_orbitals_PTR->electronegativity = 1.5;
        }
    if (atom == "W") // 6s2 5d4 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 4, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 74;
        atom_orbitals_PTR->reduced_Z = 20;
        atom_orbitals_PTR->electronegativity = 2.36;
        }
    if (atom == "Re") // 6s2 5d5 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 5, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 75;
        atom_orbitals_PTR->reduced_Z = 21;
        atom_orbitals_PTR->electronegativity = 1.9;
        }
    if (atom == "Os") // 6s2 5d6 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 6, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 76;
        atom_orbitals_PTR->reduced_Z = 22;
        atom_orbitals_PTR->electronegativity = 2.2;
        }
    if (atom == "Ir") // 6s2 5d7 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 7, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 77;
        atom_orbitals_PTR->reduced_Z = 23;
        atom_orbitals_PTR->electronegativity = 2.2;
        }
    if (atom == "Pt") // 6s2 5d8 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 9, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 78;
        atom_orbitals_PTR->reduced_Z = 24;
        atom_orbitals_PTR->electronegativity = 2.28;
        }
    if (atom == "Au") // 6s1 5d10 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 79;
        atom_orbitals_PTR->reduced_Z = 25;
        atom_orbitals_PTR->electronegativity = 2.54;
        }
    if (atom == "Hg") // 6s2 5d10 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 80;
        atom_orbitals_PTR->reduced_Z = 26;
        atom_orbitals_PTR->electronegativity = 2;
        }
    if (atom == "Tl") // 6s2 6p1 5d10 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 1, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 81;
        atom_orbitals_PTR->reduced_Z = 27;
        atom_orbitals_PTR->electronegativity = 1.62;
        }
    if (atom == "Pb") // 6s2 6p2 5d10 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 1, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 82;
        atom_orbitals_PTR->reduced_Z = 28;
        atom_orbitals_PTR->electronegativity = 2.33;
        }
    if (atom == "Bi") // 6s2 6p3 5d10 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 1, 3, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 83;
        atom_orbitals_PTR->reduced_Z = 29;
        atom_orbitals_PTR->electronegativity = 2.02;
        }
    if (atom == "Po") // 6s2 6p4 5d10 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 1, 4, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 84;
        atom_orbitals_PTR->reduced_Z = 30;
        atom_orbitals_PTR->electronegativity = 2;
        }
    if (atom == "At") // 6s2 6p5 5d10 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 1, 5, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 85;
        atom_orbitals_PTR->reduced_Z = 31;
        atom_orbitals_PTR->electronegativity = 2.2;
        }
    if (atom == "Rn") // 6s2 6p6 5d10 4f14
        {
        Quantum_numbers_to_orbitals(6, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 1, 6, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(4, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 86;
        atom_orbitals_PTR->reduced_Z = 32;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Fr") // 7s1
        {
        Quantum_numbers_to_orbitals(7, 0, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 87;
        atom_orbitals_PTR->reduced_Z = 1;
        atom_orbitals_PTR->electronegativity = 0.7;
        }
    if (atom == "Ra") // 7s2
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 88;
        atom_orbitals_PTR->reduced_Z = 2;
        atom_orbitals_PTR->electronegativity = 0.9;
        }
    if (atom == "Ac") // 7s2 6d1
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 1, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 89;
        atom_orbitals_PTR->reduced_Z = 3;
        atom_orbitals_PTR->electronegativity = 1.1;
        }
    if (atom == "Th") // 7s2 6d2
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 90;
        atom_orbitals_PTR->reduced_Z = 4;
        atom_orbitals_PTR->electronegativity = 1.3;
        }
    if (atom == "Pa") // 7s2 6d1 5f2
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 2, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 91;
        atom_orbitals_PTR->reduced_Z = 5;
        atom_orbitals_PTR->electronegativity = 1.5;
        }
    if (atom == "U") // 7s2 6d1 5f3
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 3, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 92;
        atom_orbitals_PTR->reduced_Z = 6;
        atom_orbitals_PTR->electronegativity = 1.38;
        }
    if (atom == "Np") // 7s2 6d1 5f4
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 4, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 93;
        atom_orbitals_PTR->reduced_Z = 7;
        atom_orbitals_PTR->electronegativity = 1.34;
        }
    if (atom == "Pu") // 7s2 5f6
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 6, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 94;
        atom_orbitals_PTR->reduced_Z = 8;
        atom_orbitals_PTR->electronegativity = 1.28;
        }
    if (atom == "Am") // 7s2 5f7
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 7, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 95;
        atom_orbitals_PTR->reduced_Z = 9;
        atom_orbitals_PTR->electronegativity = 1.3;
        }
    if (atom == "Cm") // 7s2 6d1 5f7
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 7, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 96;
        atom_orbitals_PTR->reduced_Z = 10;
        atom_orbitals_PTR->electronegativity = 1.3;
        }
    if (atom == "Bk") // 7s2 5f9
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 9, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 97;
        atom_orbitals_PTR->reduced_Z = 11;
        atom_orbitals_PTR->electronegativity = 1.3;
        }
    if (atom == "Cf") // 7s2 5f10
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 10, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 98;
        atom_orbitals_PTR->reduced_Z = 12;
        atom_orbitals_PTR->electronegativity = 1.3;
        }
    if (atom == "Es") // 7s2 5f11
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 11, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 99;
        atom_orbitals_PTR->reduced_Z = 13;
        atom_orbitals_PTR->electronegativity = 1.3;
        }
    if (atom == "Fm") // 7s2 5f12
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 12, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 100;
        atom_orbitals_PTR->reduced_Z = 14;
        atom_orbitals_PTR->electronegativity = 1.3;
        }
    if (atom == "Md") // 7s2 5f13
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 13, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 101;
        atom_orbitals_PTR->reduced_Z = 15;
        atom_orbitals_PTR->electronegativity = 1.3;
        }
    if (atom == "No") // 7s2 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 102;
        atom_orbitals_PTR->reduced_Z = 16;
        atom_orbitals_PTR->electronegativity = 1.3;
        }
    if (atom == "Lr") // 7s2 7p1 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(7, 1, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 103;
        atom_orbitals_PTR->reduced_Z = 17;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Rf") // 7s2 6d2 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 104;
        atom_orbitals_PTR->reduced_Z = 18;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Db") // 7s2 6d3 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 3, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 105;
        atom_orbitals_PTR->reduced_Z = 19;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Sg") // 7s2 6d4 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 4, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 106;
        atom_orbitals_PTR->reduced_Z = 20;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Bh") // 7s2 6d5 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 5, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 107;
        atom_orbitals_PTR->reduced_Z = 21;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Hs") // 7s2 6d6 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 6, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 108;
        atom_orbitals_PTR->reduced_Z = 22;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Mt") // 7s2 6d7 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 7, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 109;
        atom_orbitals_PTR->reduced_Z = 23;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Ds") // 7s2 6d8 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 8, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 110;
        atom_orbitals_PTR->reduced_Z = 24;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Rg") // 7s2 6d9 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 9, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 111;
        atom_orbitals_PTR->reduced_Z = 25;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Cn") // 7s2 6d10 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 112;
        atom_orbitals_PTR->reduced_Z = 26;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Nh") // 7s2 7p1 6d10 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(7, 1, 1, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 113;
        atom_orbitals_PTR->reduced_Z = 27;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Fl") // 7s2 7p2 6d10 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(7, 1, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 114;
        atom_orbitals_PTR->reduced_Z = 28;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Mc") // 7s2 7p3 6d10 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(7, 1, 3, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 115;
        atom_orbitals_PTR->reduced_Z = 29;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Lv") // 7s2 7p4 6d10 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(7, 1, 4, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 116;
        atom_orbitals_PTR->reduced_Z = 30;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Ts") // 7s2 7p5 6d10 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(7, 1, 5, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 117;
        atom_orbitals_PTR->reduced_Z = 31;
        atom_orbitals_PTR->electronegativity = 0;
        }
    if (atom == "Og") // 7s2 7p6 6d10 5f14
        {
        Quantum_numbers_to_orbitals(7, 0, 2, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(7, 1, 6, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(6, 2, 10, atom_orbitals_PTR);
        Quantum_numbers_to_orbitals(5, 3, 14, atom_orbitals_PTR);
        
        atom_orbitals_PTR->Z = 118;
        atom_orbitals_PTR->reduced_Z = 32;
        atom_orbitals_PTR->electronegativity = 0;
        }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Create_atomic_wavefunctions(atom_orbitals *atom_orbitals_PTR, atom_wavefunctions *atom_wavefunctions,
    unsigned int size_order,T x,T y,T z)
    { // create or modify atom_wavefunctions list of wavefunctions, probabilities and effective lenghts according to lenght multipliers
    unsigned int i, j;
    unsigned int Z;
    unsigned int count_orbitals;
    unsigned int count_electrons;
    
    count_orbitals = 0;
    count_electrons = atom_orbitals_PTR->n.size();
    Z = atom_orbitals_PTR->Z;
    
    atom_wavefunctions->lenghts.reserve(count_electrons);
    atom_wavefunctions->wavefunctions.reserve(count_electrons);
    atom_wavefunctions->probabilities.reserve(count_electrons);
    atom_wavefunctions->lenght_orders.reserve(count_electrons);
    atom_wavefunctions->x_range.reserve(count_electrons);
    atom_wavefunctions->y_range.reserve(count_electrons);
    atom_wavefunctions->z_range.reserve(count_electrons);
    atom_wavefunctions->wavefunction_lenght_multipliers.reserve(count_electrons);
    atom_wavefunctions->wavefunction_coefficients.reserve(count_electrons);
    atom_wavefunctions->effective_radius_base.reserve(count_electrons);
    atom_wavefunctions->spins.reserve(count_electrons);
    atom_wavefunctions->spin_paired.reserve(count_electrons);
    atom_wavefunctions->electron_numbers.reserve(count_electrons);
    atom_wavefunctions->bonding.reserve(count_electrons);
    atom_wavefunctions->antibonding.reserve(count_electrons);
    atom_wavefunctions->wavefunction_constraints.reserve(count_electrons);
    atom_wavefunctions->n.reserve(count_electrons);
    atom_wavefunctions->l.reserve(count_electrons);
    atom_wavefunctions->m.reserve(count_electrons);
    atom_wavefunctions->charge.reserve(count_electrons);
    atom_wavefunctions->count_electrons.reserve(count_electrons);
    atom_wavefunctions->Z.reserve(count_electrons);
    atom_wavefunctions->reduced_Z.reserve(count_electrons);
    atom_wavefunctions->x.reserve(count_electrons);
    atom_wavefunctions->y.reserve(count_electrons);
    atom_wavefunctions->z.reserve(count_electrons);
    
        
    j = 0;
    for (i = 0; i < count_electrons; i++) // Solving of spin-pairing and copying of next parameters
        {
        if (atom_orbitals_PTR->paired_with_previous[i] == false)
            j = i;
            
        if (i < count_electrons -1)
            {
            if (atom_orbitals_PTR->paired_with_previous[i] == false and atom_orbitals_PTR->paired_with_previous[i + 1] == false)
                atom_wavefunctions->spin_paired.push_back(-1);
            else
                {
                if (atom_orbitals_PTR->paired_with_previous[i] == false)
                    atom_wavefunctions->spin_paired.push_back(electron_number + 1);
                else
                    atom_wavefunctions->spin_paired.push_back(electron_number - 1);
                }
            }
        else
            {
            if (atom_orbitals_PTR->paired_with_previous[i] == false)
                atom_wavefunctions->spin_paired.push_back(-1);
            else
                atom_wavefunctions->spin_paired.push_back(electron_number - 1);
            }
        atom_wavefunctions->lenght_orders.push_back(size_order);
        atom_wavefunctions->x_range.push_back(size_order);
        atom_wavefunctions->y_range.push_back(size_order);
        atom_wavefunctions->z_range.push_back(size_order);
        atom_wavefunctions->count_electrons.push_back(count_electrons);
        atom_wavefunctions->spins.push_back(atom_orbitals_PTR->s[i]);
        atom_wavefunctions->n.push_back(atom_orbitals_PTR->n[i]);
        atom_wavefunctions->l.push_back(atom_orbitals_PTR->l[i]);
        atom_wavefunctions->m.push_back(atom_orbitals_PTR->m[i]);
        atom_wavefunctions->Z.push_back(atom_orbitals_PTR->Z);
        atom_wavefunctions->charge.push_back(atom_orbitals_PTR->charge);
        atom_wavefunctions->reduced_Z.push_back(atom_orbitals_PTR->reduced_Z);
        atom_wavefunctions->electron_numbers.push_back(electron_number);
        atom_wavefunctions->bonding.push_back(-1);
        atom_wavefunctions->antibonding.push_back(-1);
        atom_wavefunctions->pi_bonding.push_back(-1);
        atom_wavefunctions->wavefunction_constraints.push_back(0);
        atom_wavefunctions->wavefunction_coefficients.push_back(1);
        atom_wavefunctions->x.push_back(x);
        atom_wavefunctions->y.push_back(y);
        atom_wavefunctions->z.push_back(z);
        atom_wavefunctions->wavefunction_lenght_multipliers.push_back(Get_relative_Hartree_length(atom_orbitals_PTR->Z,
        atom_orbitals_PTR->n[i]));
        electron_number++;
        }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Create_bond_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions_1,
atom_wavefunctions *atom_wavefunctions_2, unsigned int count, T electronegativity_1, T electronegativity_2,
T x_difference, T y_difference, T z_difference)
    {
    unsigned int i, j, k;
    unsigned int position;
    unsigned int count_bonds; // count of the remaining electron bond pair requirements
    unsigned int count_to_process; // count of found electron pairs for bonding
    unsigned int atom_wavefunctions_1_size, atom_wavefunctions_2_size;
    unsigned int bonding_count_1, bonding_count_2;
    unsigned int unbonding_count_1, unbonding_count_2;
    unsigned int antibonding_index;
    unsigned int spin_index;
    array<unsigned int, max_electrons> indexes_1;
    array<unsigned int, max_electrons> indexes_2;
    array<unsigned int, max_electrons> indexes_3;
    array<unsigned int, max_electrons> indexes_4;
    unsigned int i_1 = 0, i_2 = 0, i_3 = 0, i_4 = 0;
    unsigned int new_l;
    
    int m;
    int preferred_m, preferred_m_begin;
    int new_m;
    bool antibonding;
    bool spin_changed;
    bool found, processed;
    
    T polarity;
    T coefficient_1, coefficient_2;
    T spins[] = {0.5, -0.5, 0.5, -0.5,  0.5, -0.5};
    
    count_bonds = count;
    
    if (abs(x_difference) >= abs(y_difference) and abs(x_difference) >= abs(z_difference)) // set preferred l
        preferred_m = 1;
    else 
        if (abs(z_difference) >= abs(y_difference))
            preferred_m = 0;
        else
            preferred_m = -1;
    
    preferred_m_begin = preferred_m;      
    atom_wavefunctions_1_size = atom_wavefunctions_1->n.size();
    atom_wavefunctions_2_size = atom_wavefunctions_2->n.size();
    
    for (k = 0; k < 3; k++) // searching electrons for bonding with various preferred m 
        {
        i_1 = 0;
        i_2 = 0;
        i_3 = 0;
        i_4 = 0;
        for (i = 0; i < atom_wavefunctions_1_size; i++) // search for posible electrons in atom 1
            {
            if (((atom_wavefunctions_1->l[i] == 0) or (atom_wavefunctions_1->m[i] > 0 and preferred_m > 0) or
            (atom_wavefunctions_1->m[i] < 0 and  preferred_m < 0) or (atom_wavefunctions_1->m[i] == 0 and  preferred_m == 0))
            and atom_wavefunctions_1->bonding[i] == -1)
                {
                if (atom_wavefunctions_1->spin_paired[i] == -1)
                    {
                    indexes_1[i_1] = i;
                    i_1++;
                    }
                else
                    {
                    indexes_3[i_3] = i;
                    i_3++;
                    }
                }
            }
        for (i = 0; i < atom_wavefunctions_2_size; i++) // search for possible electrons in atom 2
            {
            if (((atom_wavefunctions_2->l[i] == 0) or (atom_wavefunctions_2->m[i] > 0 and  preferred_m > 0) or
            (atom_wavefunctions_2->m[i] < 0 and  preferred_m < 0) or (atom_wavefunctions_2->m[i] == 0 and  preferred_m == 0))
            and atom_wavefunctions_2->bonding[i] == -1)
                {
                if (atom_wavefunctions_2->spin_paired[i] == -1)
                    {
                    indexes_2[i_2] = i;
                    i_2++;
                    }
                else
                    {
                    indexes_4[i_4] = i;
                    i_4++;
                    }
                }
            }
        if (i_1 < count_bonds) // rotating electrons by changing m if necessary and possible for atom_wavefunction_1
            for (i = 0; i < atom_wavefunctions_1_size; i++) 
                {
                if (((atom_wavefunctions_1->m[i] <= 0 and  preferred_m > 0) or
                (atom_wavefunctions_1->m[i] >= 0 and  preferred_m < 0) or (atom_wavefunctions_1->m[i] != 0 and  preferred_m == 0))
                and atom_wavefunctions_1->bonding[i] == -1)
                    { // search for posible electrons in atom 1
                    new_l = atom_wavefunctions_1->l[i];
                    for (new_m = preferred_m * new_l; ((new_m > 0 and preferred_m == 1) or (new_m == 0 and preferred_m == 0)
                    or (new_m < 0 and preferred_m == -1)); new_m -= preferred_m)
                        {
                        found = false;
                        for (j = 0; j < atom_wavefunctions_1_size; j++)
                            {
                            if (new_l == atom_wavefunctions_1->l[j] and // search for electrons with the same m
                            new_m == atom_wavefunctions_1->m[j])
                                found = true;
                            }
                        if (found == false)
                            {
                            atom_wavefunctions_1->l[i] = new_l;
                            atom_wavefunctions_1->m[i] = new_m;
                            atom_wavefunctions_1->spins[i] = 0.5;
                            indexes_1[i_1] = i;
                            i_1++;
                            break;
                            }
                        if (preferred_m == 0)
                            break;
                        }
                    }
                }
        if (i_2 < count_bonds) // rotating electrons by changing m if necessary and possible for atom_wavefunction_2
            for (i = 0; i < atom_wavefunctions_2_size; i++) 
                {
                if (((atom_wavefunctions_2->m[i] <= 0 and  preferred_m > 0) or
                (atom_wavefunctions_2->m[i] >= 0 and  preferred_m < 0) or (atom_wavefunctions_2->m[i] != 0 and  preferred_m == 0))
                and atom_wavefunctions_2->bonding[i] == -1) // search for posible electrons in atom 1
                    {
                    new_l = atom_wavefunctions_2->l[i];
                    for (new_m = preferred_m * new_l; ((new_m > 0 and preferred_m == 1) or (new_m == 0 and preferred_m == 0)
                    or (new_m < 0 and preferred_m == -1)); new_m -= preferred_m)
                        { 
                        found = false;
                        for (j = 0; j < atom_wavefunctions_2_size; j++)
                            if (new_l == atom_wavefunctions_2->l[j] and new_m == atom_wavefunctions_2->m[j])
                                found = true; // search for electrons with the same m
                                
                        if (found == false)
                            {
                            atom_wavefunctions_2->l[i] = new_l;
                            atom_wavefunctions_2->m[i] = new_m;
                            atom_wavefunctions_2->spins[i] = 0.5;
                            indexes_2[i_2] = i;
                            i_2++;
                            break;
                            }
                        if (preferred_m == 0)
                            break;
                        }
                    }
                }
        if (i_1 < count_bonds) // rotating electrons by changing m and l if necessary and possible for atom_wavefunction_1
            for (i = 0; i < atom_wavefunctions_1_size; i++) 
                {
                if (((atom_wavefunctions_1->m[i] <= 0 and  preferred_m > 0) or
                (atom_wavefunctions_1->m[i] >= 0 and  preferred_m < 0) or (atom_wavefunctions_1->m[i] != 0 and  preferred_m == 0))
                and atom_wavefunctions_1->bonding[i] == -1)
                    { // search for posible electrons in atom 1
                    for (new_l = atom_wavefunctions_1->l[i]; new_l > atom_wavefunctions_1->n[i]; new_l++)
                        for (new_m = preferred_m * new_l; ((new_m > 0 and preferred_m == 1) or (new_m == 0 and preferred_m == 0)
                        or (new_m < 0 and preferred_m == -1)); new_m -= preferred_m)
                            {
                            found = false;
                            for (j = 0; j < atom_wavefunctions_1_size; j++)
                                {
                                if (new_l == atom_wavefunctions_1->l[j] and // search for electrons with the same m
                                new_m == atom_wavefunctions_1->m[j])
                                    found = true;
                                }
                            if (found == false)
                                {
                                atom_wavefunctions_1->l[i] = new_l;
                                atom_wavefunctions_1->m[i] = new_m;
                                atom_wavefunctions_1->spins[i] = 0.5;
                                indexes_1[i_1] = i;
                                i_1++;
                                break;
                                }
                            if (preferred_m == 0)
                                break;
                            }
                    }
                }
        if (i_2 < count_bonds) // rotating electrons by changing m and l if necessary and possible for atom_wavefunction_2
            for (i = 0; i < atom_wavefunctions_2_size; i++) 
                {
                if (((atom_wavefunctions_2->m[i] <= 0 and  preferred_m > 0) or
                (atom_wavefunctions_2->m[i] >= 0 and  preferred_m < 0) or (atom_wavefunctions_2->m[i] != 0 and  preferred_m == 0))
                and atom_wavefunctions_2->bonding[i] == -1) // search for posible electrons in atom 1
                    {
                    for (new_l = atom_wavefunctions_2->l[i]; new_l > atom_wavefunctions_2->n[i]; new_l++);
                        for (new_m = preferred_m * new_l; ((new_m > 0 and preferred_m == 1) or (new_m == 0 and preferred_m == 0)
                        or (new_m < 0 and preferred_m == -1)); new_m -= preferred_m)
                            {
                            found = false;
                            for (j = 0; j < atom_wavefunctions_2_size; j++)
                                {
                                if (new_l == atom_wavefunctions_2->l[j] and new_m == atom_wavefunctions_2->m[j])
                                    found = true; // search for electrons with the same m
                                }
                            if (found == false)
                                {
                                atom_wavefunctions_2->l[i] = new_l;
                                atom_wavefunctions_2->m[i] = new_m;
                                atom_wavefunctions_2->spins[i] = 0.5;
                                indexes_2[i_2] = i;
                                i_2++;
                                break;
                                }
                            if (preferred_m == 0)
                                break;
                            }
                    }
                }
        if (count_bonds > i_1) // removal of excess electrons from orbitals for bonding for atom_wavefunction_1
            for (spin_index = 0; spin_index < 6; spin_index++)
                {
                for (i = 0; i < i_3; i++)
                    {
                    processed = false;
                    for (new_l = 0; new_l <= atom_wavefunctions_1->n[indexes_3[i]]/2; new_l++)
                        for (new_m = new_l; abs(new_m) <= new_l; new_m--)
                            { // cycles for excitation of electrons in filled orbitals
                            found = false;
                            for (j = 0; j < atom_wavefunctions_1->n.size(); j++)
                                {
                                if (atom_wavefunctions_1->l[j] == new_l and atom_wavefunctions_1->m[j] == new_m and
                                   (atom_wavefunctions_1->spins[j] == spins[spin_index] or
                                   (spin_index < 1 and atom_wavefunctions_1->spins[j] == -spins[spin_index]) or 
                                   (spin_index < 3 and atom_wavefunctions_1->bonding[j] >= 0)))
                                    found = true; 
                                } // search for 1. empty orbitals 2. half filled non-bonding orbitals 3. half-filled bonding orbitals
                            if (found == true or atom_wavefunctions_1->spin_paired[indexes_3[i]] == -1 or processed == true)
                                continue; // check avialability of orbitals and for previous processing
                        
                            for (j = 0; j < atom_wavefunctions_1->n.size(); j++) 
                                { // finding and processing electron with in the same spin-orbital
                                if (atom_wavefunctions_1->l[j] == atom_wavefunctions_1->l[indexes_3[i]] and
                                    atom_wavefunctions_1->m[j] == atom_wavefunctions_1->m[indexes_3[i]] and
                                    atom_wavefunctions_1->spins[j] != atom_wavefunctions_1->spins[indexes_3[i]])
                                    {
                                    indexes_1[i_1] = j;
                                    atom_wavefunctions_1->spins[j] = 0.5;
                                    atom_wavefunctions_1->spin_paired[j] = -1;
                                    i_1++;
                                    
                                    // Jump electron in indexes_3[i] position to new orbital
                                    atom_wavefunctions_1->l[indexes_3[i]] = new_l; 
                                    atom_wavefunctions_1->m[indexes_3[i]] = new_m;
                                    atom_wavefunctions_1->spins[indexes_3[i]] = spins[spin_index];
                                    atom_wavefunctions_1->spin_paired[indexes_3[i]] = -1;
                                    processed = true;
                                    break;
                                    }
                                }
                            
                            for (j = 0; j < atom_wavefunctions_1->n.size(); j++) 
                                { // finding and processing electron with in the orbital and test for bonding
                                if (atom_wavefunctions_1->n[j] == atom_wavefunctions_1->n[indexes_3[i]] and
                                    atom_wavefunctions_1->l[j] == atom_wavefunctions_1->l[indexes_3[i]] and
                                    atom_wavefunctions_1->m[j] == atom_wavefunctions_1->m[indexes_3[i]] and j != indexes_3[i])
                                    { // spin pairing with non-bonding electrons found in the same orbital
                                    if (atom_wavefunctions_1->bonding[k] == -1)
                                        {
                                        atom_wavefunctions_1->spin_paired[indexes_3[i]] = atom_wavefunctions_1->electron_numbers[j];
                                        atom_wavefunctions_1->spin_paired[j] = atom_wavefunctions_1->electron_numbers[indexes_3[i]];
                                        }
                                    }
                                }
                            }
                    if (count_bonds <= i_1)
                        break;
                    }
                if (count_bonds <= i_1)
                    break;
                }
        if (count_bonds > i_2) // removal of excess electrons from orbitals for bonding for atom_wavefunction_2
            for (spin_index = 0; spin_index < 6; spin_index++)
                {
                for (i = 0; i < i_4; i++)
                    {
                    processed = false;
                    for (new_l = 0; new_l <= atom_wavefunctions_2->n[indexes_4[i]]/2; new_l++)
                        for (new_m = new_l; abs(new_m) <= new_l; new_m--)
                            { // cycles for excitation of electrons in filled orbitals
                            found = false;
                            for (j = 0; j < atom_wavefunctions_2->n.size(); j++)
                                {
                                if (atom_wavefunctions_2->l[j] == new_l and atom_wavefunctions_2->m[j] == new_m and
                                   (atom_wavefunctions_2->spins[j] == spins[spin_index] or
                                   (spin_index < 1 and atom_wavefunctions_2->spins[j] == -spins[spin_index]) or 
                                   (spin_index < 3 and atom_wavefunctions_2->bonding[j] >= 0)))
                                    found = true;
                                }
                            if (found == true or atom_wavefunctions_2->spin_paired[indexes_4[i]] == -1 or processed == true)
                                continue;
                                
                            for (j = 0; j < atom_wavefunctions_2->n.size(); j++) 
                                { // finding and processing electron with in the same spin-orbital
                                if (atom_wavefunctions_2->l[j] == atom_wavefunctions_2->l[indexes_4[i]] and
                                   atom_wavefunctions_2->m[j] == atom_wavefunctions_2->m[indexes_4[i]] and
                                   atom_wavefunctions_2->spins[j] != atom_wavefunctions_2->spins[indexes_4[i]])
                                    {
                                    indexes_2[i_2] = j;
                                    atom_wavefunctions_2->spins[j] = 0.5;
                                    atom_wavefunctions_2->spin_paired[j] = -1;
                                    i_2++;
                                    
                                    // Jump electron in indexes_4[i] position to new orbital
                                    atom_wavefunctions_2->l[indexes_4[i]] = new_l; 
                                    atom_wavefunctions_2->m[indexes_4[i]] = new_m;
                                    atom_wavefunctions_2->spins[indexes_4[i]] = spins[spin_index];
                                    atom_wavefunctions_2->spin_paired[indexes_4[i]] = -1;
                                    processed = true;
                                    break;
                                    }
                                }
                            for (j = 0; j < atom_wavefunctions_2->n.size(); j++) 
                                { // finding and processing electron with in the same spin-orbital and test for bonding
                                if (atom_wavefunctions_2->n[j] == atom_wavefunctions_2->n[indexes_4[i]] and
                                    atom_wavefunctions_2->l[j] == atom_wavefunctions_2->l[indexes_4[i]] and
                                    atom_wavefunctions_2->m[j] == atom_wavefunctions_2->m[indexes_4[i]] and j != indexes_4[i])
                                    { // spin pairing with non-bonding electrons
                                    if (atom_wavefunctions_2->bonding[j] == -1)
                                        {
                                        atom_wavefunctions_2->spin_paired[indexes_4[i]] = atom_wavefunctions_2->electron_numbers[j];
                                        atom_wavefunctions_2->spin_paired[j] = atom_wavefunctions_2->electron_numbers[indexes_4[i]];
                                        }
                                    }
                                }
                            }
                    if (count_bonds <= i_2)
                        break;
                    }
                if (count_bonds <= i_2)
                    break;
                }
        if (i_1 > i_2)  // setting number of electron pairs for processing
            count_to_process = i_2;
        else
            count_to_process = i_1;
            
        if (count_to_process > count_bonds)
            count_to_process = count_bonds;
            
        for (i = 0; i < count_to_process; i++)
            { // bonding electrons are spin - paired
            atom_wavefunctions_1->spin_paired[indexes_1[i]] = atom_wavefunctions_2->electron_numbers[indexes_2[i]];
            atom_wavefunctions_2->spin_paired[indexes_2[i]] = atom_wavefunctions_1->electron_numbers[indexes_1[i]];
            antibonding = false;
            spin_changed = false;
            if (atom_wavefunctions_1->spins[indexes_1[i]] == atom_wavefunctions_2->spins[indexes_2[i]])
                {
                atom_wavefunctions_2->spins[indexes_2[i]] = - atom_wavefunctions_2->spins[indexes_2[i]]; // changing spin if necessary
                spin_changed = true;
                }
            for (j = 0; j < atom_wavefunctions_2->n.size(); j++) // unpaired electron to antibonding orbital
                if (atom_wavefunctions_2->n[j] == atom_wavefunctions_2->n[indexes_2[i]] and 
                atom_wavefunctions_2->l[j] == atom_wavefunctions_2->l[indexes_2[i]] and
                atom_wavefunctions_2->m[j] == atom_wavefunctions_2->m[indexes_2[i]] and
                j != indexes_2[i])
                    {
                    if (atom_wavefunctions_2->bonding[j] == -1)
                        {
                        if (spin_changed == true) // spin swap
                            {
                            atom_wavefunctions_2->spins[j] = - atom_wavefunctions_2->spins[j];
                            spin_changed == false;
                            }
                        atom_wavefunctions_2->spin_paired[j] = -1;
                        atom_wavefunctions_2->antibonding[j] = atom_wavefunctions_1->electron_numbers[indexes_1[i]];
                        antibonding = true;
                        antibonding_index = j;
                        }
                    }
            for (j = 0; j < atom_wavefunctions_1->n.size(); j++) // next unpaired electron to antibonding orbital
                if (atom_wavefunctions_1->n[j] == atom_wavefunctions_1->n[indexes_1[i]] and
                atom_wavefunctions_1->l[j] == atom_wavefunctions_1->l[indexes_1[i]] and
                atom_wavefunctions_1->m[j] == atom_wavefunctions_1->m[indexes_1[i]] and
                j != indexes_1[i])
                    if (atom_wavefunctions_1->bonding[j] == -1)
                        {
                        if (antibonding == true)
                            {
                            atom_wavefunctions_1->spin_paired[j] = atom_wavefunctions_2->electron_numbers[antibonding_index];
                            atom_wavefunctions_1->antibonding[j] = atom_wavefunctions_2->electron_numbers[indexes_2[i]];
                            atom_wavefunctions_2->spin_paired[antibonding_index] = atom_wavefunctions_1->electron_numbers[j];
                            }
                        else
                            {
                            atom_wavefunctions_1->spin_paired[j] = -1;
                            atom_wavefunctions_1->antibonding[j] = -1;
                            }
                        }
            atom_wavefunctions_1->bonding[indexes_1[i]] = atom_wavefunctions_2->electron_numbers[indexes_2[i]];
            atom_wavefunctions_2->bonding[indexes_2[i]] = atom_wavefunctions_1->electron_numbers[indexes_1[i]];
            if (preferred_m_begin != preferred_m and atom_wavefunctions_1->l[indexes_1[i]] != 0
            and atom_wavefunctions_2->l[indexes_2[i]] != 0)
                {
                atom_wavefunctions_1->pi_bonding[indexes_1[i]] = atom_wavefunctions_2->electron_numbers[indexes_2[i]];
                atom_wavefunctions_2->pi_bonding[indexes_2[i]] = atom_wavefunctions_1->electron_numbers[indexes_1[i]];
                }
            count_bonds--;
            polarity = 1 - exp(-(electronegativity_1 - electronegativity_2)*(electronegativity_1 - electronegativity_2)/4);
            // Linus Pauling formula
            if ((sqrt(x_difference * x_difference + y_difference * y_difference + z_difference * z_difference)
            /(atom_wavefunctions_1->n[indexes_1[i]] + atom_wavefunctions_2->n[indexes_2[i]])) > 1) // cut-off function
                {
                polarity = polarity * exp(-(sqrt(x_difference * x_difference + y_difference * y_difference + z_difference * z_difference)
                /(atom_wavefunctions_1->n[indexes_1[i]] + atom_wavefunctions_2->n[indexes_2[i]])) + 1);
                }
            if (electronegativity_1 > electronegativity_2)
                {
                coefficient_1 = sqrt(1 + polarity);
                coefficient_2 = sqrt(1 - polarity);
                }
            else
                {
                coefficient_1 = sqrt(1 - polarity);
                coefficient_2 = sqrt(1 + polarity);
                }
            atom_wavefunctions_1->wavefunction_coefficients[indexes_1[i]] = coefficient_1;
            atom_wavefunctions_2->wavefunction_coefficients[indexes_2[i]] = coefficient_2;
            }
        switch (preferred_m)
            {
            case -1:
                preferred_m = 0;
                break;
            case 0:
                preferred_m = 1;
                break;
            case 1:
                preferred_m = -1;
                break;
            }
        }
    if (count_bonds > 0)
        return(-1); // Found not sufficient amount of electrons for all bonds
        
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Sum_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions_1, atom_wavefunctions *atom_wavefunctions_2)
    {
    unsigned int i;
    unsigned int size_1, size_2;
    int charge, charge_2;
    
    T x;
    T y;
    T z;
    T Z, Z_2; // for adding inner valence shell by noble gas atoms
    T reduced_Z, reduced_Z_2;
    T count_electrons, count_electrons_2;
    Z_2 = 0;
    reduced_Z_2 = 0;
    count_electrons_2 = 0;
    size_1 = atom_wavefunctions_1->n.size();
    size_2 = atom_wavefunctions_2->n.size();
    
    if (size_1 > 0 and size_2 > 0)
        {
        x = atom_wavefunctions_2->x[0];
        y = atom_wavefunctions_2->y[0];
        z = atom_wavefunctions_2->z[0];
        Z = atom_wavefunctions_2->Z[0];
        reduced_Z = atom_wavefunctions_2->reduced_Z[0];
        count_electrons = atom_wavefunctions_2->count_electrons[0];
        charge = atom_wavefunctions_2->charge[0];
        Z_2 = Z;
        charge_2 = charge;
        reduced_Z_2 = reduced_Z;
        count_electrons_2 = count_electrons;
        
        for (i = size_1 - 1; i >= 0 ; i--)
            {
            if (atom_wavefunctions_1->x[i] == x and atom_wavefunctions_1->y[i] == y and atom_wavefunctions_1->z[i] == z)
                {
                if (Z > atom_wavefunctions_1->Z[i])
                    atom_wavefunctions_1->Z[i] = Z;
                // Backward update of reduced_Z and count_eleectrons
                atom_wavefunctions_1->reduced_Z[i] = atom_wavefunctions_1->reduced_Z[i] + reduced_Z; 
                atom_wavefunctions_1->count_electrons[i] = atom_wavefunctions_1->count_electrons[i] + count_electrons;
                atom_wavefunctions_1->charge[i] = atom_wavefunctions_1->charge[i] + charge;
                
                Z_2 = atom_wavefunctions_1->Z[i];
                reduced_Z_2 = atom_wavefunctions_1->reduced_Z[i];
                count_electrons_2 = atom_wavefunctions_1->count_electrons[i];
                charge_2 = atom_wavefunctions_1->charge[i];
                }
            if (i == 0 or atom_wavefunctions_1->x[i] != x or atom_wavefunctions_1->y[i] != y or atom_wavefunctions_1->z[i] != z)
                break;
            }
        }
    for (i = 0; i < size_2; i++)
        {
        atom_wavefunctions_1->lenght_orders.push_back(atom_wavefunctions_2->lenght_orders[i]);
        atom_wavefunctions_1->x_range.push_back(atom_wavefunctions_2->x_range[i]);
        atom_wavefunctions_1->y_range.push_back(atom_wavefunctions_2->y_range[i]);
        atom_wavefunctions_1->z_range.push_back(atom_wavefunctions_2->z_range[i]);
        atom_wavefunctions_1->wavefunction_lenght_multipliers.push_back(atom_wavefunctions_2->wavefunction_lenght_multipliers[i]);
        atom_wavefunctions_1->wavefunction_coefficients.push_back(atom_wavefunctions_2->wavefunction_coefficients[i]);
        atom_wavefunctions_1->spins.push_back(atom_wavefunctions_2->spins[i]);
        atom_wavefunctions_1->spin_paired.push_back(atom_wavefunctions_2->spin_paired[i]);
        atom_wavefunctions_1->bonding.push_back(atom_wavefunctions_2->bonding[i]);
        atom_wavefunctions_1->antibonding.push_back(atom_wavefunctions_2->antibonding[i]);
        atom_wavefunctions_1->pi_bonding.push_back(atom_wavefunctions_2->pi_bonding[i]);
        atom_wavefunctions_1->wavefunction_constraints.push_back(atom_wavefunctions_2->wavefunction_constraints[i]);
        atom_wavefunctions_1->n.push_back(atom_wavefunctions_2->n[i]);
        atom_wavefunctions_1->l.push_back(atom_wavefunctions_2->l[i]);
        atom_wavefunctions_1->m.push_back(atom_wavefunctions_2->m[i]);
        if ((atom_wavefunctions_2->x[i] == x and atom_wavefunctions_2->y[i] == y and atom_wavefunctions_2->z[i] == z) and size_1 > 0)
            {
            atom_wavefunctions_1->Z.push_back(Z_2); // for adding inner valence shell by noble gas atoms
            atom_wavefunctions_1->reduced_Z.push_back(reduced_Z_2);
            atom_wavefunctions_1->count_electrons.push_back(count_electrons_2);
            atom_wavefunctions_1->charge.push_back(charge_2);
            }
        else
            {
            atom_wavefunctions_1->Z.push_back(atom_wavefunctions_2->Z[i]);
            atom_wavefunctions_1->reduced_Z.push_back(atom_wavefunctions_2->reduced_Z[i]);
            atom_wavefunctions_1->count_electrons.push_back(atom_wavefunctions_2->count_electrons[i]);
            atom_wavefunctions_1->charge.push_back(atom_wavefunctions_2->charge[i]);
            }
        atom_wavefunctions_1->electron_numbers.push_back(atom_wavefunctions_2->electron_numbers[i]);
        atom_wavefunctions_1->x.push_back(atom_wavefunctions_2->x[i]);
        atom_wavefunctions_1->y.push_back(atom_wavefunctions_2->y[i]);
        atom_wavefunctions_1->z.push_back(atom_wavefunctions_2->z[i]);
        }
    return(0);
    }
// End of section 4 - generating list of electrons, section 5 - generating matrices of integrals and Fock matrices
template <typename T>
int Slater_basis_set_calculations<T>::Create_nuclear_atraction_integral_matrix(T* matrix, T* nucleuses,
unsigned int order, atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j, k;
    unsigned int count_orbitals;
    unsigned int sum_electrons;
    unsigned int ind_nuc_size;
    array<unsigned int, max_electrons> index; // Index of electrons positions forcomputing nuclear atraction integrals
    array<unsigned int, max_electrons + 1> ind_nuc; // Index of nucleuses
    
    array<T*, max_electrons> probabilities;
    array<unsigned int, max_electrons> lenght_orders;
    array<unsigned int, max_electrons> x_range;
    array<unsigned int, max_electrons> y_range;
    array<unsigned int, max_electrons> z_range;
    array<T, max_electrons> wavefunction_coefficients;
    array<T, max_electrons> wavefunction_lenght_multipliers;
    array<T, max_electrons> effective_radius_base;
    array<unsigned int, max_electrons> count_electrons;
    array<unsigned int, max_electrons> Z;
    array<T, max_electrons> x;
    array<T, max_electrons> y;
    array<T, max_electrons> z;
    array<unsigned int, max_electrons> n;
    array<int, max_electrons> bonding;
    array<int, max_electrons> spin_paired;
    array<T, max_electrons> spins;
    array<unsigned int, max_electrons> wavefunction_constraints;
    
    bool restriction = true;
    
    sum_electrons = atom_wavefunctions->n.size();
    count_orbitals = 0;
    
    if (sum_electrons > max_electrons)
        return(-1);
    
    copy(atom_wavefunctions->probabilities.begin(), atom_wavefunctions->probabilities.end(), probabilities.begin());
    copy(atom_wavefunctions->lenght_orders.begin(), atom_wavefunctions->lenght_orders.end(), lenght_orders.begin());
    copy(atom_wavefunctions->x_range.begin(), atom_wavefunctions->x_range.end(), x_range.begin());
    copy(atom_wavefunctions->y_range.begin(), atom_wavefunctions->y_range.end(), y_range.begin());
    copy(atom_wavefunctions->z_range.begin(), atom_wavefunctions->z_range.end(), z_range.begin());
    copy(atom_wavefunctions->wavefunction_coefficients.begin(), atom_wavefunctions->wavefunction_coefficients.end(),
    wavefunction_coefficients.begin());
    copy(atom_wavefunctions->wavefunction_lenght_multipliers.begin(),
    atom_wavefunctions->wavefunction_lenght_multipliers.end(), wavefunction_lenght_multipliers.begin());
    copy(atom_wavefunctions->effective_radius_base.begin(), atom_wavefunctions->effective_radius_base.end(),
    effective_radius_base.begin());
    copy(atom_wavefunctions->count_electrons.begin(), atom_wavefunctions->count_electrons.end(),
    count_electrons.begin());
    copy(atom_wavefunctions->reduced_Z.begin(), atom_wavefunctions->reduced_Z.end(), Z.begin());
    copy(atom_wavefunctions->x.begin(), atom_wavefunctions->x.end(), x.begin());
    copy(atom_wavefunctions->y.begin(), atom_wavefunctions->y.end(), y.begin());
    copy(atom_wavefunctions->z.begin(), atom_wavefunctions->z.end(), z.begin());
    copy(atom_wavefunctions->n.begin(), atom_wavefunctions->n.end(), n.begin());
    copy(atom_wavefunctions->spins.begin(), atom_wavefunctions->spins.end(), spins.begin());
    copy(atom_wavefunctions->spin_paired.begin(), atom_wavefunctions->spin_paired.end(), spin_paired.begin());
    copy(atom_wavefunctions->bonding.begin(), atom_wavefunctions->bonding.end(), bonding.begin());
    copy(atom_wavefunctions->wavefunction_constraints.begin(), atom_wavefunctions->wavefunction_constraints.end(),
    wavefunction_constraints.begin());
    
    // filling matrices with a zeros, values of frozen orbitals are remained
    for (i = 0; i < order; i++)
        if (wavefunction_constraints[i] == 0 or iterations == 0)
            for (j = 0; j < order; j++)
                {
                matrix[i * order + j] = 0;
                nucleuses[i * order + j] = 0;
                }
    // closed-shell basis set method optimalization code
    for (i = 0; i < sum_electrons ; i++) // for restricted and unrestricted basis set method
        {
        if (spin_paired[i] == -1)
            {
            restriction = false;
            }
        }
    for (i = 0; i < sum_electrons ; i++)
        if (wavefunction_coefficients[i] != 0)
            {
            if (i == 0 and (wavefunction_constraints[i] == 0 or iterations == 0))
                {
                index[count_orbitals] = i;
                count_orbitals++;
                }
            if (i > 0)
                {
                if ((spins[i] == 0.5 or (bonding[i] >= 0 or restriction == false))
                and (wavefunction_constraints[i] == 0 or iterations == 0))
                    {
                    index[count_orbitals] = i;
                    count_orbitals++;
                    }
                }
            }
    i = 0;
    j = 1;
    ind_nuc[0] = 0;
    for (i = 1; i < sum_electrons; i++)
        if (x[i] != x[i-1] or y[i] != y[i-1] or z[i] != z[i-1])
            {
            ind_nuc[j] = i;
            j++;
            }
    ind_nuc_size = j;
    ind_nuc[ind_nuc_size] = order;
    
    // multithreading code
    for (i = 0; i < ind_nuc_size; i++)
        {
        # pragma omp parallel
            {
            # pragma omp for
            for (j = 0; j < count_orbitals; j++)
                {
                if (((x[index[j]] - x[ind_nuc[i]] != 0) or (y[index[j]] - y[ind_nuc[i]] != 0) or
                (z[index[j]] - z[ind_nuc[i]] != 0)) and ((n[index[j]]/(Z[index[j]] *
                wavefunction_lenght_multipliers[index[j]])) > (vector_lenght/lenght_orders[index[j]])))
                    {
                    Integrate_Integral_nucleus_atraction(probabilities[index[j]], matrix + (index[j] * (1 + order)),
                    nucleuses_distances + (index[j] * order + ind_nuc[i]),lenght_orders[index[j]],
                    x_range[index[j]], y_range[index[j]], z_range[index[j]], x[index[j]] - x[ind_nuc[i]],
                    y[index[j]] - y[ind_nuc[i]], z[index[j]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
                    }
                }
            }
        for (j = 0; j < count_orbitals; j++)
            {
            if (not(((x[index[j]] - x[ind_nuc[i]] != 0) or (y[index[j]] - y[ind_nuc[i]] != 0) or
            (z[index[j]] - z[ind_nuc[i]] != 0)) and ((n[index[j]]/(Z[index[j]] *
            wavefunction_lenght_multipliers[index[j]])) > (vector_lenght/lenght_orders[index[j]]))))
                {
                Integral_nucleus_atraction(effective_radius_base[index[j]],
                wavefunction_lenght_multipliers[index[j]], matrix + (index[j] * (order + 1)),
                nucleuses_distances + (index[j] * order + ind_nuc[i]), Z[ind_nuc[i]]);
                }
            }
        for(j = 0; j < order; j++) // copying values into colums of nucleuses
            {
            nucleuses[i + (j * order)] = matrix[j * (order + 1)] * wavefunction_coefficients[j]
            * wavefunction_coefficients[j]; // copy diagonal of matrix
            
            for (k = 0; k < i; k++)
                nucleuses[i + (j * order)] = nucleuses[i + (j * order)] - nucleuses[k + (j * order)];
                
            } // subtract previous values in the row
        }
    // End of multithreading code
    // closed-shell basis set method optimalization code
    if (restriction == true) // copying values in restricted basis_set_method
        for (i = 0; i < order; i++)
            {
            if (spins[i] == -0.5 and bonding[i] == -1 and spin_paired[i] >= 0)
                {
                for (j = 0; j < ind_nuc_size; j++)
                    nucleuses_distances[i * order + ind_nuc[j]] = nucleuses_distances[spin_paired[i] * order + ind_nuc[j]];
            
                matrix[i * (1 + order)] = matrix[spin_paired[i] * (1 + order)];
                
                for (j = 0; j < order; j++)
                    nucleuses[i * order + j] = nucleuses[spin_paired[i] * order + j];
                }
            }
    // end closed-shell basis set method optimalization code
    for (i = 0; i < order; i++) // including the wavefunctions coefficients for linear combination
        {
        matrix[i * (1 + order)] = matrix[i * (1 + order)] * wavefunction_coefficients[i] * wavefunction_coefficients[i];
        
        for (j = 0; j < ind_nuc_size; j++)
            for (k = ind_nuc[j] + 1; k < ind_nuc[j + 1]; k++)
                nucleuses_distances[i * order + k] = nucleuses_distances[i * order + ind_nuc[j]];
        }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Create_coulombic_integral_matrix(T* matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions,
small_atom_wavefunctions *small_atom_wavefunctions)
    {
    unsigned int i, j, k;
    unsigned int count_electrons, count_orbitals;
    unsigned int index_size, index_2_size;
    unsigned int small_atom_wavefunctions_size;
    array<unsigned int, max_electrons> index;
    array<unsigned int, max_electrons> index_2_array;
    
    array<T, max_electrons> wavefunction_coefficients;
    array<T, max_electrons> effective_radius_base;
    array<T, max_electrons> wavefunction_lenght_multipliers;
    array<unsigned int, max_electrons> lenght_orders;
    array<unsigned int, max_electrons> x_range;
    array<unsigned int, max_electrons> y_range;
    array<unsigned int, max_electrons> z_range;
    array<int, max_electrons> bonding;
    array<unsigned int, max_electrons> n;
    array<unsigned int, max_electrons> l;
    array<T, max_electrons> spins;
    array<int, max_electrons> spin_paired;
    array<T, max_electrons> x;
    array<T, max_electrons> y;
    array<T, max_electrons> z;
    array<T*, max_electrons> small_probabilities;
    array<unsigned int, max_electrons> wavefunction_constraints;
    T distance, radius_1, radius_2;
    bool restriction;
    bool spin_bonded;
    
    count_electrons = atom_wavefunctions->n.size();
    small_atom_wavefunctions_size = small_atom_wavefunctions->n.size();
    count_orbitals = 0;
    restriction = true;
    
    copy(atom_wavefunctions->wavefunction_coefficients.begin(), atom_wavefunctions->wavefunction_coefficients.end(),
    wavefunction_coefficients.begin());
    copy(atom_wavefunctions->effective_radius_base.begin(), atom_wavefunctions->effective_radius_base.end(),
    effective_radius_base.begin());
    copy(atom_wavefunctions->wavefunction_lenght_multipliers.begin(),
    atom_wavefunctions->wavefunction_lenght_multipliers.end(), wavefunction_lenght_multipliers.begin());
    copy(small_atom_wavefunctions->lenght_orders.begin(), small_atom_wavefunctions->lenght_orders.end(),
    lenght_orders.begin());
    copy(atom_wavefunctions->x_range.begin(), atom_wavefunctions->x_range.end(), x_range.begin());
    copy(atom_wavefunctions->y_range.begin(), atom_wavefunctions->y_range.end(), y_range.begin());
    copy(atom_wavefunctions->z_range.begin(), atom_wavefunctions->z_range.end(), z_range.begin());
    copy(atom_wavefunctions->bonding.begin(), atom_wavefunctions->bonding.end(), bonding.begin());
    copy(atom_wavefunctions->n.begin(), atom_wavefunctions->n.end(), n.begin());
    copy(atom_wavefunctions->l.begin(), atom_wavefunctions->l.end(), l.begin());
    copy(atom_wavefunctions->spins.begin(), atom_wavefunctions->spins.end(), spins.begin());
    copy(atom_wavefunctions->spin_paired.begin(), atom_wavefunctions->spin_paired.end(), spin_paired.begin());
    copy(atom_wavefunctions->x.begin(), atom_wavefunctions->x.end(), x.begin());
    copy(atom_wavefunctions->y.begin(), atom_wavefunctions->y.end(), y.begin());
    copy(atom_wavefunctions->z.begin(), atom_wavefunctions->z.end(), z.begin());
    copy(small_atom_wavefunctions->probabilities.begin(), small_atom_wavefunctions->probabilities.end(),
    small_probabilities.begin());
    copy(atom_wavefunctions->wavefunction_constraints.begin(), atom_wavefunctions->wavefunction_constraints.end(),
    wavefunction_constraints.begin());
    
    // filling matrices with a zeros, values of frozen orbitals are remained
    for (i = 0; i < order; i++)
        for (j = 0; j < order; j++)
            if (wavefunction_constraints[i] == 0 or wavefunction_constraints[j] == 0 or iterations == 0)
                matrix[i * order + j] = 0;
    
    for (i = 0; i < count_electrons; i++) // Using regression curve for s1 - s1 integrals
        {
        for (j = i + 1; j < count_electrons; j++)
            if ((n[i] == 1 and n[j] == 1) or
            ((n[i] == 1 or n[j] == 1) and (x[j] - x[i] == 0 and y[j] - y[i] == 0 and z[j] - z[i] == 0))
            or (spin_paired[i] == j and (x[j] - x[i] == 0 and y[j] - y[i] == 0 and z[j] - z[i] == 0)))
                {
                radius_1 = effective_radius_base[i] * wavefunction_lenght_multipliers[i];
                radius_2 = effective_radius_base[j] * wavefunction_lenght_multipliers[j];
                distance = sqrt(((x[j] - x[i]) * (x[j] - x[i])) + ((y[j] - y[i]) * (y[j] - y[i])) + ((z[j] - z[i])
                * (z[j] - z[i]))) * Hartree_lenght;
                if (spin_paired[i] == j and Helium_correlation_energy == true)
                    spin_bonded = true;
                else
                    spin_bonded = false;
                Integral_coulombic(radius_1, radius_2, distance, matrix + (i + (j * order)), spin_bonded);
                }
        }
    // closed-shell basis set method optimalization code
    for (i = 0; i < count_electrons ; i++)
        {
        if (atom_wavefunctions->spin_paired[i] == -1)// for molecules an open-shell systems
            restriction = false;
        }
    for (i = 0; i < count_electrons; i++) // list of electrons for calculation
        if (wavefunction_coefficients[i] != 0)
            {
            if (restriction == false)
                {
                index[count_orbitals] = i;
                count_orbitals++;
                }
            else
                {
                if ((bonding[i] >= 0) or spin_paired[i] > i)
                    { // for restricted basis set method only in bonding orbitals both electrons
                    index[count_orbitals] = i;
                    count_orbitals++;
                    }
                }
            }
    index_size = count_orbitals;
    for (i = 0; i < count_orbitals; i++)
        {
        k = 0;
        for (j = i + 1; j < index_size; j++)
            {
            if (matrix[index[i] + (index[j] * order)] != 0)
                continue;
            
            index_2_array[k] = index[j];
            k++;
            }
        if (restriction == true and bonding[index[i]] == -1)
            {
            if (matrix[index[i] + (spin_paired[index[i]] * order)] != 0)
                continue;
            
            index_2_array[k] = spin_paired[index[i]];
            k++;
            }
        index_2_size = k;
        // multithreading code
        #pragma omp parallel
            {
            #pragma omp for
            for (j = 0; j < index_2_size; j++)
                {
                Integrate_Integral_coulombic(small_probabilities[index[i]],
                small_probabilities[index_2_array[j]], &matrix[index[i] + (index_2_array[j] * order)],
                lenght_orders[index_2_array[j]], x[index_2_array[j]] - x[index[i]], y[index_2_array[j]] - y[index[i]],
                z[index_2_array[j]] - z[index[i]], x_range[index[i]], x_range[index_2_array[j]],
                y_range[index[i]], y_range[index_2_array[j]], z_range[index[i]], z_range[index_2_array[j]]);
                }
            }
        } // End of multithreading code
    // closed-shell basis set method optimalization code 
    if (restriction == true) // copy calculated values for electrons in pairs with same wavefunction pairs
        for (i = 0; i < index_size; i++)
            for (j = i + 1; j < index_size; j++)
                { // prevention rewriting within electron pairs
                if (matrix[spin_paired[index[i]] + (index[j] * order)] == 0 and bonding[index[i]] == -1) 
                    matrix[spin_paired[index[i]] + (index[j] * order)] = matrix[index[i] + (index[j] * order)];
                                    
                if (matrix[index[i]              + ((spin_paired[index[j]]) * order)] == 0 and bonding[index[j]] == -1)
                    matrix[index[i]              + ((spin_paired[index[j]]) * order)] = matrix[index[i] + (index[j] * order)];
                  
                if (matrix[spin_paired[index[i]] + ((spin_paired[index[j]]) * order)] == 0 and bonding[index[i]] == -1
                and bonding[index[j]] == -1)
                    matrix[spin_paired[index[i]] + ((spin_paired[index[j]]) * order)] = matrix[index[i] + (index[j] * order)];
                }
    // end of closed-shell basis set method optimalization code 
    for (i = 0; i < order; i++) // including the wavefunctions coefficients for linear combination
        for (j = 0; j < i; j++)
            {
            if (bonding[i] != j)
                matrix[j + (i * order)] = matrix[j + (i * order)] * wavefunction_coefficients[i] *
                wavefunction_coefficients[i] * wavefunction_coefficients[j] * wavefunction_coefficients[j];
            else
                { // for orbitals with wavefunction_coefficients > 1 are coulombic and exchange energies cancelled
                if (wavefunction_coefficients[j] > 1)
                    matrix[j + (i * order)] = matrix[j + (i * order)] * wavefunction_coefficients[i] * wavefunction_coefficients[i];
                if (wavefunction_coefficients[i] > 1)
                    matrix[j + (i * order)] = matrix[j + (i * order)] * wavefunction_coefficients[j] * wavefunction_coefficients[j];
                }
            }
    for (i = 0; i < order; i++) // copy calculated effective_lenghts upper diagonal
        for (j = 0; j < i; j++)
            matrix[(j * order) + i] = matrix[(i * order) + j];
            
    for (i = 0; i < order; i++) // copying 0 to diagonal
        matrix[i * (order + 1)] = 0;
        
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Create_overlap_integral_matrix(T* matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j, k;
    unsigned int count_electrons;
    unsigned int count_orbitals;
    unsigned int index_size;
    unsigned int index_2_size;
    
    array<unsigned int, max_electrons> index;
    array<unsigned int, max_electrons> index_2_array;
    array<T*, max_electrons> wavefunctions;
    array<unsigned int, max_electrons> lenght_orders;
    array<unsigned int, max_electrons> x_range;
    array<unsigned int, max_electrons> y_range;
    array<unsigned int, max_electrons> z_range;
    array<T, max_electrons> wavefunction_coefficients;
    array<T, max_electrons> x;
    array<T, max_electrons> y;
    array<T, max_electrons> z;
    array<T, max_electrons> spins;
    array<unsigned int, max_electrons> l;
    array<int, max_electrons> m;
    array<int, max_electrons> spin_paired;
    array<int, max_electrons> bonding;
    array<unsigned int, max_electrons> wavefunction_constraints;
    bool restriction;
    
    restriction = true;
    count_electrons = atom_wavefunctions->n.size();
    count_orbitals = 0;
    
    copy(atom_wavefunctions->wavefunctions.begin(), atom_wavefunctions->wavefunctions.end(),
    wavefunctions.begin());
    copy(atom_wavefunctions->lenght_orders.begin(), atom_wavefunctions->lenght_orders.end(),
    lenght_orders.begin());
    copy(atom_wavefunctions->x_range.begin(), atom_wavefunctions->x_range.end(), x_range.begin());
    copy(atom_wavefunctions->y_range.begin(), atom_wavefunctions->y_range.end(), y_range.begin());
    copy(atom_wavefunctions->z_range.begin(), atom_wavefunctions->z_range.end(), z_range.begin());
    copy(atom_wavefunctions->wavefunction_coefficients.begin(), atom_wavefunctions->wavefunction_coefficients.end(),
    wavefunction_coefficients.begin());
    copy(atom_wavefunctions->x.begin(), atom_wavefunctions->x.end(), x.begin());
    copy(atom_wavefunctions->y.begin(), atom_wavefunctions->y.end(), y.begin());
    copy(atom_wavefunctions->z.begin(), atom_wavefunctions->z.end(), z.begin());
    copy(atom_wavefunctions->spins.begin(), atom_wavefunctions->spins.end(), spins.begin());
    copy(atom_wavefunctions->l.begin(), atom_wavefunctions->l.end(), l.begin());
    copy(atom_wavefunctions->m.begin(), atom_wavefunctions->m.end(), m.begin());
    copy(atom_wavefunctions->spin_paired.begin(), atom_wavefunctions->spin_paired.end(), spin_paired.begin());
    copy(atom_wavefunctions->bonding.begin(), atom_wavefunctions->bonding.end(), bonding.begin());
    copy(atom_wavefunctions->wavefunction_constraints.begin(), atom_wavefunctions->wavefunction_constraints.end(),
    wavefunction_constraints.begin());
    
    // filling matrices with a zeros, values of frozen orbitals are remained
    for (i = 0; i < order; i++)
        for (j = 0; j < order; j++)
            if (wavefunction_constraints[i] == 0 or wavefunction_constraints[j] == 0 or iterations == 0)
                matrix[i * order + j] = 0;
    
    // closed-shell basis set method optimalization code
    for (i = 0; i < count_electrons ; i++)
        {
        if (atom_wavefunctions->spin_paired[i] == -1)// for molecules an open-shell systems
            restriction = false;
        }
    for (i = 0; i < count_electrons; i++) // list of electrons for calculation
        if (wavefunction_coefficients[i] != 0)
            {
            if (restriction == false)
                {
                index[count_orbitals] = i;
                count_orbitals++;
                }
            else
                {
                if ((bonding[i] >= 0) or spin_paired[i] > i)
                    { // for restricted basis set method only the bonding orbitals both electrons
                    index[count_orbitals] = i;
                    count_orbitals++;
                    }
                }
            }
    index_size = count_orbitals;
    for (i = 0; i < count_orbitals; i++)
        {
        k = 0;
        if (restriction == false)
            for (j = i + 1; j < index_size; j++)
                {
                if (x[index[i]] == x[index[j]] and y[index[i]] == y[index[j]] and z[index[i]] == z[index[j]] and
                   (l[index[i]] != l[index[j]] or m[index[i]] != m[index[j]]))
                    continue; // Electrons with different l or m in the same atom have zero overlap integral.
                    
                if ((x[index[i]] != x[index[j]] or y[index[i]] != y[index[j]] or z[index[i]] != z[index[j]]) and 
                   ((m[index[i]] != m[index[j]] and l[index[i]] == l[index[j]]) or l[index[i]] == 0 or l[index[j]] == 0) and
                   (  (((l[index[i]] == 1 and m[index[i]] == 1) or (l[index[j]] == 1 and m[index[j]] == 1)) and
                   x[index[i]] - x[index[j]] == 0)
                   or (((l[index[i]] == 1 and m[index[i]] == 0) or (l[index[j]] == 1 and m[index[j]] == 0)) and
                   z[index[i]] - z[index[j]] == 0)
                   or (((l[index[i]] == 1 and m[index[i]] == -1) or (l[index[j]] == 1 and m[index[j]] == -1)) and
                   y[index[i]] - y[index[j]] == 0)))
                    continue; // p, d and f electrons in different atoms with zero overlap integral
                
                if (wavefunction_constraints[i] != 0 and wavefunction_constraints[j] != 0 and iterations > 0)
                    continue;
                    
                index_2_array[k] = index[j];
                k++;
                }
        else
            {
            for (j = i + 1; j < index_size; j++)
                if (bonding[index[j]] >= 0)
                    {
                    if (x[index[i]] == x[index[j]] and y[index[i]] == y[index[j]] and z[index[i]] == z[index[j]] and
                       (l[index[i]] != l[index[j]] or m[index[i]] != m[index[j]]))
                        continue; // Electrons with different l or m in the same atom have zero overlap integral.
                        
                    if ((x[index[i]] != x[index[j]] or y[index[i]] != y[index[j]] or z[index[i]] != z[index[j]]) and 
                    ((m[index[i]] != m[index[j]] and l[index[i]] == l[index[j]]) or l[index[i]] == 0 or l[index[j]] == 0) and
                    (  (((l[index[i]] == 1 and m[index[i]] == 1) or (l[index[j]] == 1 and m[index[j]] == 1)) and
                    x[index[i]] - x[index[j]] == 0)
                    or (((l[index[i]] == 1 and m[index[i]] == 0) or (l[index[j]] == 1 and m[index[j]] == 0)) and
                    z[index[i]] - z[index[j]] == 0)
                    or (((l[index[i]] == 1 and m[index[i]] == -1) or (l[index[j]] == 1 and m[index[j]] == -1)) and
                    y[index[i]] - y[index[j]] == 0)))
                        continue; // p electrons in different atoms with zero overlap integral
                    
                    index_2_array[k] = index[j];
                    k++;
                    }
            }
        index_2_size = k;
        // end of closed-shell basis set method optimalization code    
        // multithreading code
        #pragma omp parallel
            {
            #pragma omp for
            for (j = 0; j < index_2_size; j++)
                {
                Integral_overlap(wavefunctions[index[i]], wavefunctions[index_2_array[j]],
                &matrix[index[i] + (index_2_array[j] * order)], lenght_orders[index_2_array[j]],
                x_range[index[i]], x_range[index_2_array[j]], y_range[index[i]], y_range[index_2_array[j]],
                z_range[index[i]], z_range[index_2_array[j]],
                x[index_2_array[j]] - x[index[i]], y[index_2_array[j]] - y[index[i]], z[index_2_array[j]] - z[index[i]]);
                }
            }
        } // end of multithreading code
    // closed-shell basis set method optimalization code 
    if (restriction == true) // copy calculated values for electrons in pairs with same wavefunction pairs
        for (i = 0; i < index_size; i++)
            for (j = i + 1; j < index_size; j++)
                { // prevention rewriting within electron pairs
                if (matrix[spin_paired[index[i]] + (index[j] * order)] == 0 and bonding[index[i]] == -1) 
                    matrix[spin_paired[index[i]] + (index[j] * order)] = matrix[index[i] + (index[j] * order)];
                    
                if (matrix[index[i]              + ((spin_paired[index[j]]) * order)] == 0 and bonding[index[j]] == -1)
                    matrix[index[i]              + ((spin_paired[index[j]]) * order)] = matrix[index[i] + (index[j] * order)];
                          
                if (matrix[spin_paired[index[i]] + ((spin_paired[index[j]]) * order)] == 0 and bonding[index[i]] == -1
                and bonding[index[j]] == -1)
                    matrix[spin_paired[index[i]] + ((spin_paired[index[j]]) * order)] = matrix[index[i] + (index[j] * order)];
                }
    // end of closed-shell basis set method optimalization code 
    for (i = 0; i < order; i++) // copy 1 absolute value for electrons in pairs with same wavefunction pairs
        if ((restriction == true) and (spin_paired[i] >= i) and  (bonding[i] == -1))
        matrix[i + spin_paired[i] * order] = 1;
              
    for (i = 0; i < order; i++) // including wavefunction coefficients of linear combination and spins
        for (j = i + 1; j < order; j++)
            matrix[(j * order) + i] = matrix[(j * order) + i] * wavefunction_coefficients[i] * wavefunction_coefficients[j]
            * spins[i] * spins[j] * 4;
            
    for (i = 0; i < order; i++) // copy calculated effective_lenghts upper diagonal
        for (j = 0; j < i; j++)
            matrix[(j * order) + i] = matrix[(i * order) + j];
            
    for (i = 0; i < order; i++) // copying 1 to diagonal
        matrix[i * (order + 1)] = 1;

    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Calculate_resonance_integral_matrix(T* overlap_matrix, T* overlap_effective_lenght_integral_matrix,
T* resonance_integral_matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions,
small_atom_wavefunctions *small_atom_wavefunctions)
    {
    unsigned int i, j, k;
    unsigned int count_electrons, count_small_wavefunctions;
    unsigned int count_interactions;
    
    array<T, max_electrons> spins;
    array<T, max_electrons> x;
    array<T, max_electrons> y;
    array<T, max_electrons> z;
    array<unsigned int, max_electrons> n;
    array<unsigned int, max_electrons> electron_numbers;
    array<T, max_electrons> effective_radius_base;
    array<T, max_electrons> wavefunction_lenght_multipliers;
    array<T, max_electrons> wavefunction_coefficients;
    array<int, max_electrons> bonding;
    array<int, max_electrons> antibonding;
    array<unsigned int, max_electrons> small_electron_numbers;
    array<unsigned int, max_electrons> small_lenght_orders;
    array<T*, max_electrons> small_wavefunctions;
    array<T, max_electrons> small_x;
    array<T, max_electrons> small_y;
    array<T, max_electrons> small_z;
    array<unsigned int, max_electrons> index;
    array<unsigned int, max_electrons> wavefunction_constraints;
    unsigned int index_size;
    T radius, radius_1, radius_2;
    T effective_lenght;
    T constant;
    T overlap;
    T pow_overrlaps;
    
    count_electrons = atom_wavefunctions->n.size();
    count_small_wavefunctions = small_atom_wavefunctions->electron_numbers.size();
    
    
    copy(atom_wavefunctions->spins.begin(), atom_wavefunctions->spins.end(), spins.begin());
    copy(atom_wavefunctions->x.begin(), atom_wavefunctions->x.end(), x.begin());
    copy(atom_wavefunctions->y.begin(), atom_wavefunctions->y.end(), y.begin());
    copy(atom_wavefunctions->z.begin(), atom_wavefunctions->z.end(), z.begin());
    copy(atom_wavefunctions->n.begin(), atom_wavefunctions->n.end(), n.begin());
    
    copy(atom_wavefunctions->electron_numbers.begin(), atom_wavefunctions->electron_numbers.end(),
    electron_numbers.begin());
    copy(atom_wavefunctions->effective_radius_base.begin(), atom_wavefunctions->effective_radius_base.end(),
    effective_radius_base.begin());
    copy(atom_wavefunctions->wavefunction_lenght_multipliers.begin(),
    atom_wavefunctions->wavefunction_lenght_multipliers.end(), wavefunction_lenght_multipliers.begin());
    copy(atom_wavefunctions->wavefunction_coefficients.begin(), atom_wavefunctions->wavefunction_coefficients.end(),
    wavefunction_coefficients.begin());
    copy(atom_wavefunctions->bonding.begin(), atom_wavefunctions->bonding.end(), bonding.begin());
    copy(atom_wavefunctions->antibonding.begin(), atom_wavefunctions->antibonding.end(), antibonding.begin());
    
    copy(small_atom_wavefunctions->electron_numbers.begin(), small_atom_wavefunctions->electron_numbers.end(),
    small_electron_numbers.begin());
    copy(small_atom_wavefunctions->lenght_orders.begin(), small_atom_wavefunctions->lenght_orders.end(),
    small_lenght_orders.begin());
    copy(small_atom_wavefunctions->wavefunctions.begin(), small_atom_wavefunctions->wavefunctions.end(),
    small_wavefunctions.begin());
    copy(small_atom_wavefunctions->x.begin(), small_atom_wavefunctions->x.end(), small_x.begin());
    copy(small_atom_wavefunctions->y.begin(), small_atom_wavefunctions->y.end(), small_y.begin());
    copy(small_atom_wavefunctions->z.begin(), small_atom_wavefunctions->z.end(), small_z.begin());
    copy(atom_wavefunctions->wavefunction_constraints.begin(), atom_wavefunctions->wavefunction_constraints.end(),
    wavefunction_constraints.begin());
    
    constant = e*e/(4*Pi*E0);
    
    // filling matrices with a zeros, values of frozen orbitals are remained
    for (i = 0; i < order; i++)
        for (j = 0; j < order; j++)
            if (wavefunction_constraints[i] == 0 or wavefunction_constraints[j] == 0 or iterations == 0)
                {
                overlap_effective_lenght_integral_matrix[i * order + j] = 0;
                resonance_integral_matrix[i * order + j] = 0;
                }
    
    for (i = 0; i < order; i++)
        {
        k = 0;
        for (j = 0; j < i; j++)
            {
            if ((overlap_matrix[(i * order) + j] != 0) and (bonding[i] >= 0 or antibonding[i] >= 0) and
            (bonding[j] >= 0 or antibonding[j] >= 0) and (x[i] != x[j] or y[i] != y[j] or z[i] != z[j]))
                { // Combination of all bonding orbitals replace hybridization
                if (n[i] == 1 and n[j] == 1)
                    { // From overlap lenght integration and linear regression curve
                    radius_1 = effective_radius_base[i]/wavefunction_lenght_multipliers[i];
                    radius_2 = effective_radius_base[j]/wavefunction_lenght_multipliers[j];
                    if (radius_1 >= radius_2)
                        radius = radius_1 + (Phi + 1)/3.00 * radius_2;
                    else
                        radius = radius_2 + (Phi + 1)/3.00 * radius_1;  
                    effective_lenght = radius * (1.481178577 - (0.700616171 * abs(overlap_matrix[(i * order) + j]))
                    + (0.224068020 * abs(overlap_matrix[(i * order) + j]) * abs(overlap_matrix[(i * order) + j])));
                    overlap_effective_lenght_integral_matrix[(i * order) +j] = effective_lenght;
                    }
                else
                    { // Integrating of effective lenght of wavefunction overlap
                    // for resonance integral from small wavefunction set
                    if (wavefunction_constraints[i] == 0 or wavefunction_constraints[j] == 0 or iterations == 0)
                        {
                        index[k] = j;
                        k++;
                        }
                    }
                }
            }
        index_size = k;
        // multithreading code
        #pragma omp parallel
            {
            #pragma omp for
            for (j = 0; j < index_size; j++)
                {
                Integrate_Integral_overlap(small_wavefunctions[i], small_wavefunctions[index[j]],
                &overlap_effective_lenght_integral_matrix[(i * order) + index[j]], small_lenght_orders[i],
                small_x[index[j]] - small_x[i], small_y[index[j]] - small_y[i], small_z[index[j]] - small_z[i]);
                }
            }
        } // end of multithreading code
    for (i = 0; i < order; i++) // copy calculated effective_lenghts upper diagonal
        for (j = 0; j < i; j++)
            overlap_effective_lenght_integral_matrix[(j * order) + i] = overlap_effective_lenght_integral_matrix[(i * order) + j];
        
    for (i = 0; i < order; i++)
        for (j = 0; j < i; j++)
            if ((overlap_matrix[(i * order) + j] != 0) and (bonding[i] >= 0 or antibonding[i] >= 0) and
            (bonding[j] >= 0 or antibonding[j] >= 0) and (x[i] != x[j] or y[i] != y[j] or z[i] != z[j]))
                {
                resonance_integral_matrix[(i * order) + j] = constant/overlap_effective_lenght_integral_matrix[(i * order) + j] *
                abs(overlap_matrix[(i * order) + j] * overlap_matrix[(i * order) + j]) * spins[i] * spins[j] * 4;
                } // correction for hybridization of orbitals 
    
    for (i = 0; i < order; i++) // copy calculated values upper diagonal
        for (j = 0; j < i; j++)
            resonance_integral_matrix[(j * order) + i] = resonance_integral_matrix[(i * order) + j];
    
    for (i = 0; i < order; i++) // correction for hybridization of orbitals
        {
        pow_overrlaps = 0;
        for (j = 0; j < order; j++)
            if (resonance_integral_matrix[(i * order) + j] != 0)
                pow_overrlaps = pow_overrlaps + (overlap_matrix[(i * order) + j] * overlap_matrix[(i * order) + j]);
        
        for (j = 0; j < order; j++)
            if (resonance_integral_matrix[(i * order) + j] != 0)
                {
                resonance_integral_matrix[(i * order) + j] = resonance_integral_matrix[(i * order) + j] *
                (overlap_matrix[(i * order) + j] * overlap_matrix[(i * order) + j])/pow_overrlaps;
                }
        }
    for (i = 0; i < order; i++) // correction for  localizing bonds - row-wise
        {
        count_interactions = 0;
        for (j = 0; j < order; j++)
            if (resonance_integral_matrix[(i * order) + j] != 0)
                count_interactions++;
        
        if (count_interactions == 1)
            for (j = 0; j < order; j++)
                if (resonance_integral_matrix[(i * order) + j] != 0)
                   resonance_integral_matrix[(i * order) + j] = resonance_integral_matrix[(i * order) + j];
        }
    for (i = 0; i < order; i++) // correction for multiple delocalized bonds - column-wise
        {
        count_interactions = 0;
        for (j = 0; j < order; j++)
            if (resonance_integral_matrix[(j * order) + i] != 0)
                count_interactions++;
        
        if (count_interactions == 1)
            for (j = 0; j < order; j++)
                if (resonance_integral_matrix[(j * order) + i] != 0)
                    resonance_integral_matrix[(j * order) + i] = resonance_integral_matrix[(j * order) + i];
        }
    for (i = 0; i < order; i++) // correction for polar bonds
        for (j = 0; j < i; j++)
            {
            if (wavefunction_coefficients[i] < 1)
                resonance_integral_matrix[(i * order) + j] = resonance_integral_matrix[(i * order) + j] * wavefunction_coefficients[i];
            
            if (wavefunction_coefficients[j] < 1)
                resonance_integral_matrix[(i * order) + j] = resonance_integral_matrix[(i * order) + j] * wavefunction_coefficients[j];
            }
    for (i = 0; i < order; i++) // second copy calculated values upper diagonal
        for (j = 0; j < i; j++)
            resonance_integral_matrix[(j * order) + i] = resonance_integral_matrix[(i * order) + j];
    
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Calculate_kinetic_integral_matrix(T* matrix, unsigned int order,
atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j, k;
    unsigned int count_electrons, count_orbitals;
    unsigned int index_size, index_2_size;
    array<unsigned int, max_electrons> index;
    array<unsigned int, max_electrons> index_2_array;
    array<int, max_electrons> spin_paired;
    
    array<T*, max_electrons> Gradients;
    array<T, max_electrons> wavefunction_coefficients;
    array<T, max_electrons> effective_radius_base;
    array<T, max_electrons> wavefunction_lenght_multipliers;
    array<unsigned int, max_electrons> lenght_orders;
    array<unsigned int, max_electrons> x_range;
    array<unsigned int, max_electrons> y_range;
    array<unsigned int, max_electrons> z_range;
    array<int, max_electrons> bonding;
    array<unsigned int, max_electrons> n;
    array<unsigned int, max_electrons> l;
    array<T, max_electrons> spins;
    array<T, max_electrons> x;
    array<T, max_electrons> y;
    array<T, max_electrons> z;
    array<T, max_electrons> Z;
    array<unsigned int, max_electrons> wavefunction_constraints; /*
    T pixel_lenght;
    
    pixel_lenght = T(vector_lenght)/T(order) * Hartree_lenght;
    
    T distance, radius_1, radius_2; */
    bool restriction;
    count_electrons = atom_wavefunctions->n.size();
    count_orbitals = 0;
    restriction = true;
    
    copy(atom_wavefunctions->Gradients.begin(), atom_wavefunctions->Gradients.end(),
    Gradients.begin());
    copy(atom_wavefunctions->wavefunction_coefficients.begin(), atom_wavefunctions->wavefunction_coefficients.end(),
    wavefunction_coefficients.begin());
    copy(atom_wavefunctions->effective_radius_base.begin(), atom_wavefunctions->effective_radius_base.end(),
    effective_radius_base.begin());
    copy(atom_wavefunctions->wavefunction_lenght_multipliers.begin(),
    atom_wavefunctions->wavefunction_lenght_multipliers.end(), wavefunction_lenght_multipliers.begin());
    copy(atom_wavefunctions->lenght_orders.begin(), atom_wavefunctions->lenght_orders.end(), lenght_orders.begin());
    copy(atom_wavefunctions->x_range.begin(), atom_wavefunctions->x_range.end(), x_range.begin());
    copy(atom_wavefunctions->y_range.begin(), atom_wavefunctions->y_range.end(), y_range.begin());
    copy(atom_wavefunctions->z_range.begin(), atom_wavefunctions->z_range.end(), z_range.begin());
    copy(atom_wavefunctions->bonding.begin(), atom_wavefunctions->bonding.end(), bonding.begin());
    copy(atom_wavefunctions->n.begin(), atom_wavefunctions->n.end(), n.begin());
    copy(atom_wavefunctions->l.begin(), atom_wavefunctions->l.end(), l.begin());
    copy(atom_wavefunctions->spins.begin(), atom_wavefunctions->spins.end(), spins.begin());
    copy(atom_wavefunctions->spin_paired.begin(), atom_wavefunctions->spin_paired.end(), spin_paired.begin());
    copy(atom_wavefunctions->x.begin(), atom_wavefunctions->x.end(), x.begin());
    copy(atom_wavefunctions->y.begin(), atom_wavefunctions->y.end(), y.begin());
    copy(atom_wavefunctions->z.begin(), atom_wavefunctions->z.end(), z.begin());
    copy(atom_wavefunctions->Z.begin(), atom_wavefunctions->Z.end(), Z.begin());
    copy(atom_wavefunctions->wavefunction_constraints.begin(), atom_wavefunctions->wavefunction_constraints.end(),
    wavefunction_constraints.begin());
    
    // filling matrices with a zeros, values of frozen orbitals are remained
    for (i = 0; i < order; i++)
        for (j = 0; j < order; j++)
            if (wavefunction_constraints[i] == 0 or wavefunction_constraints[j] == 0 or iterations == 0)
                matrix[i * order + j] = 0;

    for (i = 0; i < count_electrons; i++) // Fill the diagonal one-electron kinetic energies
        {
        matrix[i * (1 + order)] = (h * h)/(8 * Pi * Pi * me * Hartree_lenght * Hartree_lenght) *
        (Z[i] * Z[i] * wavefunction_lenght_multipliers[i])/(n[i] * n[i]) * wavefunction_coefficients[i] 
        * wavefunction_coefficients[i];
        }
    // closed-shell basis set method optimalization code
    for (i = 0; i < count_electrons ; i++)
        {
        if (atom_wavefunctions->spin_paired[i] == -1)// for molecules an open-shell systems
            restriction = false;
        }
    for (i = 0; i < count_electrons; i++) // list of electrons for calculation
        if (wavefunction_coefficients[i] != 0)
            {
            if (restriction == false)
                {
                index[count_orbitals] = i;
                count_orbitals++;
                }
            else
                {
                if ((bonding[i] >= 0) or spin_paired[i] > i)
                    { // for restricted basis set method only in bonding orbitals both electrons
                    index[count_orbitals] = i;
                    count_orbitals++;
                    }
                }
            }
    index_size = count_orbitals;
    for (i = 0; i < count_orbitals; i++)
        {
        k = 0;
        for (j = i + 1; j < index_size; j++)
            {
            if (matrix[index[i] + (index[j] * order)] != 0)
                continue;
            
            index_2_array[k] = index[j];
            k++;
            }
        if (restriction == true and bonding[index[i]] == -1)
            {
            if (matrix[index[i] + (spin_paired[index[i]] * order)] != 0)
                continue;
            
            index_2_array[k] = spin_paired[index[i]];
            k++;
            }
        index_2_size = k;
        // multithreading code
        #pragma omp parallel
            {
            #pragma omp for
            for (j = 0; j < index_2_size; j++)
                {
                Integral_kinetic(Gradients[index[i]],
                Gradients[index_2_array[j]], &matrix[index[i] + (index_2_array[j] * order)],
                lenght_orders[index_2_array[j]], x[index_2_array[j]] - x[index[i]], y[index_2_array[j]] - y[index[i]],
                z[index_2_array[j]] - z[index[i]], x_range[index[i]], x_range[index_2_array[j]],
                y_range[index[i]], y_range[index_2_array[j]], z_range[index[i]], z_range[index_2_array[j]]);
                }
            }
        } // End of multithreading code
    // closed-shell basis set method optimalization code 
    if (restriction == true) // copy calculated values for electrons in pairs with same wavefunction pairs
        for (i = 0; i < index_size; i++)
            for (j = i + 1; j < index_size; j++)
                { // prevention rewriting within electron pairs
                if (matrix[spin_paired[index[i]] + (index[j] * order)] == 0 and bonding[index[i]] == -1) 
                    matrix[spin_paired[index[i]] + (index[j] * order)] = matrix[index[i] + (index[j] * order)];
                                    
                if (matrix[index[i]              + ((spin_paired[index[j]]) * order)] == 0 and bonding[index[j]] == -1)
                    matrix[index[i]              + ((spin_paired[index[j]]) * order)] = matrix[index[i] + (index[j] * order)];
                  
                if (matrix[spin_paired[index[i]] + ((spin_paired[index[j]]) * order)] == 0 and bonding[index[i]] == -1
                and bonding[index[j]] == -1)
                    matrix[spin_paired[index[i]] + ((spin_paired[index[j]]) * order)] = matrix[index[i] + (index[j] * order)];
                }
    // end of closed-shell basis set method optimalization code 
    for (i = 0; i < order; i++) // including the wavefunctions coefficients for linear combination
        for (j = 0; j < i; j++)
            { // Correct two electron kinetic energie according to lenghts
            matrix[j + (i * order)] = matrix[j + (i * order)] * wavefunction_coefficients[i] * wavefunction_coefficients[i] *
            wavefunction_coefficients[j] * wavefunction_coefficients[j];
            matrix[i + (j * order)] = matrix[j + (i * order)]; // Copy values upper diagonal
            }
    for (i = 0; i < order; i++) // copy calculated effective_lenghts upper diagonal
        for (j = 0; j < i; j++)
            matrix[(j * order) + i] = matrix[(i * order) + j];
        
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Calculate_basis_set_matrix(T* nuclear_atraction_integral_matrix, T* coulombic_integral_matrix,
T* resonance_integral_matrix, T* kinetic_integral_matrix, T* basis_set_matrix, 
unsigned int order, atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j;
    
    T spin_orbit_energy;
    T B;
    T distance_nucleuses;
    T radius, radius_1, radius_2;
    T effective_lenght;
    array<T, max_electrons> effective_radius_base;
    array<T, max_electrons> wavefunction_lenght_multipliers;
    array<T, max_electrons> potential_energy;
    
    array<unsigned int, max_electrons> Z;
    array<unsigned int, max_electrons> x;
    array<unsigned int, max_electrons> y;
    array<unsigned int, max_electrons> z;
    array<unsigned int, max_electrons> n;
    array<unsigned int, max_electrons> l;
    array<int, max_electrons> m;
    array<T, max_electrons> spins;
    array<int, max_electrons> spin_paired;
    
    for (i = 0; i < order; i++)
        {
        effective_radius_base[i] = atom_wavefunctions->effective_radius_base[i];
        wavefunction_lenght_multipliers[i] = atom_wavefunctions->wavefunction_lenght_multipliers[i];
        Z[i] = atom_wavefunctions->Z[i];
        x[i] = atom_wavefunctions->x[i];
        y[i] = atom_wavefunctions->y[i];
        z[i] = atom_wavefunctions->z[i];
        n[i] = atom_wavefunctions->n[i];
        l[i] = atom_wavefunctions->l[i];
        m[i] = atom_wavefunctions->m[i];
        spins[i] = atom_wavefunctions->spins[i];
        spin_paired[i] = atom_wavefunctions->spin_paired[i];
        potential_energy[i] = 0;
        }
    for (i = 0; i < order * order; i++)
        basis_set_matrix[i] = coulombic_integral_matrix[i] + resonance_integral_matrix[i] + kinetic_integral_matrix[i];
        
    for (i = 0; i < order; i++)
        basis_set_matrix[i * (1 + order)] = basis_set_matrix[i * (1 + order)] + nuclear_atraction_integral_matrix[i * (1 + order)];
    // adding the spin-orbit coupling
    for (i = 0; i < order; i++)
        {
        if (spin_paired[i] == -1)
            {
            for (j = 0; j < order; j++)
                potential_energy[i] = kinetic_integral_matrix[(i * order) + j];
                
            for (j = 0; j < order; j++)
                {
                if (l[j] != 0 and i != j)
                    {
                    distance_nucleuses = sqrt((x[i] - x[j] * x[i] - x[j]) + (y[i] - y[j] * y[i] - y[j]) + (z[i] - z[j] * z[i] - z[j]));
                    distance_nucleuses = abs(distance_nucleuses * Hartree_lenght);
                    
                    radius_1 = effective_radius_base[i]/wavefunction_lenght_multipliers[i];
                    radius_2 = effective_radius_base[j]/wavefunction_lenght_multipliers[j];
                    if (radius_1 >= radius_2)
                        radius = radius_1 + (Phi - 1) * radius_2;
                    else
                        radius = radius_2 + (Phi - 1) * radius_1;
    
                    effective_lenght = sqrt((radius * radius) + (distance_nucleuses * distance_nucleuses));
                       
                    B = Orbital_magnetic_field(potential_energy[i], effective_lenght, l[j]);
                    spin_orbit_energy = Spin_moment_energy(spins[i], B);
                    basis_set_matrix[(i * order) + j] = basis_set_matrix[(i * order) + j] + spin_orbit_energy;
                    basis_set_matrix[(j * order) + i] = basis_set_matrix[(j * order) + i] + spin_orbit_energy;
                    }
                }
            }
        }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Calculate_corr_basis_set_matrix(T* basis_set_matrix, T* correction_matrix,
T* corr_basis_set_matrix, unsigned int order)
    {
    unsigned int i;
    
    for (i = 0; i < order * order; i++)
        corr_basis_set_matrix[i] = basis_set_matrix[i] + correction_matrix[i];
    
    return(0);
    }
template <typename T>
T Slater_basis_set_calculations<T>::Solve_basis_set_matrix(T* basis_set_matrix, T* overlap_integral_matrix, unsigned int order,
vector<T>* values, atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j, k;
    array<unsigned int, max_electrons> indexes_2;
    
    T alpha, beta, E;
    T multiplier_constant, multiplier;
    array<T, max_electrons> Hamiltonian_parts;
    array<T, max_electrons> beta_sum_rows;
    array<T, max_electrons> overlaps;
    array<T, max_electrons> Rows;
    T Hamiltonian = 0;
    T new_old_iteration_ratio[] =  {0.50, 0.50,
    0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250,
    0.166, 0.166, 0.166, 0.166, 0.166, 0.166, 0.166, 0.166,
    0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
    0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100,
    0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833,
    0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833, 0.0833,
    0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714,
    0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714, 0.0714,};
    
    unsigned int* n = atom_wavefunctions->n.data();
    unsigned int* l = atom_wavefunctions->l.data();
    int* m = atom_wavefunctions->m.data();
    T* spins = atom_wavefunctions->spins.data();
    T* x = atom_wavefunctions->x.data();
    T* y = atom_wavefunctions->y.data();
    T* z = atom_wavefunctions->z.data();
    int* bonding = atom_wavefunctions->bonding.data();
    int* spin_paired = atom_wavefunctions->spin_paired.data();
    T* wavefunction_lenght_multipliers = atom_wavefunctions->wavefunction_lenght_multipliers.data();
    unsigned int* constraints = atom_wavefunctions->wavefunction_constraints.data();
    unsigned int* Z = atom_wavefunctions->Z.data();
    unsigned int* reduced_Z = atom_wavefunctions->reduced_Z.data();
    unsigned int* count_electrons = atom_wavefunctions->count_electrons.data();
    int* charge = atom_wavefunctions->charge.data();
    bool restriction;
    
    alpha = 0; // for diagonal elements
    beta = 0; // for non-diagonal elements
    values->clear();
    values->reserve(2 * order);
    electron_spectra.clear();
    electron_spectra.reserve(max_electrons * 2);
    restriction = true;
    
    for (i = 0; i < order; i++)
        {
        overlaps[i] = 0;
        Hamiltonian_parts[i] = 0;
        beta_sum_rows[i] = 0;
        if (atom_wavefunctions->spin_paired[i] == -1)
            restriction = false;
        }
    for (i = 0; i < order; i++)
        {
        for (j = i; j < order; j++)
            Hamiltonian_parts[i] = Hamiltonian_parts[i] + basis_set_matrix[(i * order) + j];
        }
    for (i = 0; i < order; i++) // finding alpha value
        if (abs(basis_set_matrix[i * (order + 1)]) > abs(alpha))
            alpha = basis_set_matrix[i * (order + 1)];
    
    for (i = 0; i < order; i++) // finding beta value
        for (j = i + 1; j < order; j++)
            if (abs(basis_set_matrix[i + (j * order)]) > abs(beta))
                beta = basis_set_matrix[i + (j * order)];
                
    for (i = 0; i < order; i++) // diagonal transition to shape (alpha - alpha max) / beta
        {
        basis_set_matrix[i * (order + 1)] = basis_set_matrix[i * (order + 1)] - alpha;
        if (beta != 0) 
            basis_set_matrix[i * (order + 1)] = basis_set_matrix[i * (order + 1)]/beta;
        }
    if (beta != 0) // non-diagonal transition to shape beta / beta max
        for (i = 0; i < order; i++)
            for (j = 0; j < order; j++)
                if (i != j)
                    basis_set_matrix[i + (j * order)] = basis_set_matrix[i + (j * order)]/beta;
    
    for (i = 0; i < order; i++)
        {
        for (j = 0; j < order; j++)
            beta_sum_rows[i] = beta_sum_rows[i] + basis_set_matrix[(i * order) + j];
        }
    
    basis_set_Determinant_solver(order, basis_set_matrix); // main calculation of determinants array
    for (i = 0; i < order; i++)
        {
        if (determinants.size() > 0)
            E = alpha + (determinants[0] * beta) + (beta_sum_rows[i] * beta);
            
        else
            E = alpha + beta_sum_rows[i] * beta;
            
        Rows[i] = E; // determinants into Rows sorted  according basis_set matrix rows
        }
    for (i = 0; i < order; i++) // calculating Hamiltonian
        if ((not (isnan(Hamiltonian_parts[i])) and (not isinf(Hamiltonian_parts[i])))) // Check for NaN and inf values
            Hamiltonian = Hamiltonian + Hamiltonian_parts[i];
        
    for (i = 0; i < order; i++) // copying to vectors
        {// algoritm of generating new wavefunctions cause doubling of values
        if (atom_wavefunctions->bonding[i] >= 0)
            {
            values->push_back(-abs((Rows[i] + Rows[bonding[i]])/2));
            electron_spectra.push_back(abs((Rows[i] + Rows[bonding[i]])/2));
            }
        else
            {
            values->push_back(Rows[i]);
            electron_spectra.push_back(-Rows[i]);
            }
        }
    for (i = 0; i < order; i++) // variational generation of new wavefunction lenght coefficients for next iteration from roots of matrix
        {
        if ((Z[i] - charge[i]))
            multiplier_constant = 1.00/sqrt(Z[i] - charge[i]);
            
        if (Rows[i] != 0 and -correction_matrix[i * (1 + order)] != Rydberg_energy(Z[i], n[i]) and constraints[i] == 0)
            {
            multiplier = pow(wavefunction_lenght_multipliers[i], 1 - new_old_iteration_ratio[Z[i] - charge[i]]) *
            pow(sqrt(abs(Rydberg_energy(Z[i], n[i]) /Rows[i] * multiplier_constant)), new_old_iteration_ratio[Z[i] - charge[i]])
            * Get_relative_Hartree_length(Z[i], n[i]);
            if (multiplier > 1.00/1024 and multiplier < 1024)
                wavefunction_lenght_multipliers[i] = multiplier;
            }
        }
    sort(electron_spectra.begin(),electron_spectra.end());
    return(Hamiltonian);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Generate_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions,
small_atom_wavefunctions *small_atom_wavefunctions, unsigned int size_order, bool allocate, bool compute_densities)
    { // create or modify atom_wavefunctions list of wavefunctions, probabilities and effective lenghts according to lenght multipliers
    unsigned int i, j, k;
    unsigned int count_orbitals;
    unsigned int count_electrons;
    unsigned int wavefunction_size;
    unsigned int spin_paired_electron_index;
    array<unsigned int, max_electrons> n_i;
    array<unsigned int, max_electrons> l_i;
    array<int, max_electrons> m_i;
    array<unsigned int, max_electrons> Z;
    array<unsigned int, max_electrons> Z_i;
    array<unsigned int, max_electrons> x_range;
    array<unsigned int, max_electrons> y_range;
    array<unsigned int, max_electrons> z_range;
    
    array<T*, max_electrons> pointers_to_wavefunctions;
    array<T*, max_electrons> pointers_to_probabilities;
    array<T*, max_electrons> pointers_to_Gradients;
    T* pointer_to_lenghts;
    T* pointer_to_wavefunction;
    T* pointer_to_probability;
    T* pointer_to_Gradient;
    T* small_lenghts;
    T* small_relative_lenghts;
    array<T*, max_electrons> small_wavefunctions;
    array<T*, max_electrons> small_probabilities;
    unsigned int small_lenght_order;
    unsigned int small_wavefunction_size;
    unsigned int small_atom_wavefunctions_size;
    array<unsigned int, max_electrons>  small_electron_numbers;
    array<unsigned int, max_electrons>  small_lenght_orders;
    array<unsigned int, max_electrons>  small_x_range;
    array<unsigned int, max_electrons>  small_y_range;
    array<unsigned int, max_electrons>  small_z_range;
    
    T multiplier;
    array<unsigned int, max_electrons + 1> index;
    array<T, max_electrons> multipliers;
    array<T, max_electrons> multiplier_array;
    array<T, max_electrons> effective_radius_base_array;
    array<T, max_electrons> spins;
    array<int, max_electrons> spin_paired;
    array<int, max_electrons> bonding;
    array<unsigned int, max_electrons> constraints;
    array<unsigned int, max_electrons> electron_numbers;
    
    bool restriction;
    bool non_s1_system;
    
    count_orbitals = 0;
    count_electrons = atom_wavefunctions->n.size();
    
    if (count_electrons > max_electrons)
        return(-1);
        
    wavefunction_size = ((2 * size_order) + 1) * ((2 * size_order) + 1) * ((2 * size_order) + 1);
    small_lenght_order = (sqrt(size_order) * small_wavefunction_size_factor);
    small_wavefunction_size = (2 * small_lenght_order + 1) * (2 * small_lenght_order + 1) * (2 * small_lenght_order + 1);
    
    non_s1_system = false;
    restriction = true;
    
    for (i = 0; i < count_electrons; i++) // Switch between restricted and unrestricted basis set method
        {
        if (atom_wavefunctions->spin_paired[i] == -1)
            restriction = false;
        
        Z[i] = atom_wavefunctions->Z[i];
        spins[i] = atom_wavefunctions->spins[i];
        spin_paired[i] = atom_wavefunctions->spin_paired[i];
        bonding[i] = atom_wavefunctions->bonding[i];
        constraints[i] = atom_wavefunctions->wavefunction_constraints[i];
        electron_numbers[i] = atom_wavefunctions->electron_numbers[i];
        x_range[i] = atom_wavefunctions->x_range[i];
        y_range[i] = atom_wavefunctions->y_range[i];
        z_range[i] = atom_wavefunctions->z_range[i];
        }
    // closed-shell basis set method optimalization code
    j = 0;
    if (restriction == true)
        {
        for (i = 0; i < count_electrons; i++)
            {
            if ((constraints[i] == 0) and (spins[i] == 0.5 or bonding[i] >= 0))
                {
                index[j] = i;
                j++;
                }
            }
        }
    else
        {
        for (i = 0; i < count_electrons; i++)
            {
            index[j] = i;
            j++;
            }
        }
    count_orbitals = j;
    index[count_orbitals] = count_electrons;
    for (i = 0; i < count_orbitals; i++)
        {
        n_i[i] = atom_wavefunctions->n[index[i]];
        l_i[i] = atom_wavefunctions->l[index[i]];
        m_i[i] = atom_wavefunctions->m[index[i]];
        Z_i[i] = atom_wavefunctions->Z[index[i]];
        if (n_i[i] > 1)
            non_s1_system = true;
        }
    if (allocate == true) // generating lists of lenghts for wavefunctions, wavefunctions, probabilities and lenghts multipliers
        {
        for (i = 0; i < count_electrons; i++)
            {
            pointers_to_wavefunctions[i] = nullptr;
            pointers_to_probabilities[i] = nullptr;
            pointers_to_Gradients[i] = nullptr;
            small_wavefunctions[i] = nullptr;
            small_probabilities[i] = nullptr;
            }
        
        pointer_to_lenghts = nullptr;
        small_lenghts = nullptr;
        small_relative_lenghts = nullptr;
        try
            {
            pointer_to_lenghts = new T[wavefunction_size];
            atom_wavefunctions->lenghts.push_back(pointer_to_lenghts);
            if (non_s1_system == true)
                {
                small_lenghts = new T[small_wavefunction_size];
                small_relative_lenghts = new T[small_wavefunction_size * small_wavefunction_size];
                small_atom_wavefunctions->lenghts.push_back(small_lenghts);
                small_atom_wavefunctions->relative_lenghts.push_back(small_relative_lenghts);
                }
            }
        catch (int)
            {
            if (pointer_to_lenghts != nullptr)
                delete[] pointer_to_lenghts;
            
            if (small_lenghts != nullptr)
                delete[] small_lenghts;
            
            if (small_relative_lenghts != nullptr)
                delete[] small_relative_lenghts;
            
            return(-1);
            }
        // Generating lenghts 3D cubes for wavefunction calculations
        Wavefunction_lenghts_generate(pointer_to_lenghts, size_order);
        // Generating small_lenghts 3D cubes for wavefunction calculations
        if (non_s1_system == true or s1_memory_optimization == false) {
            Wavefunction_lenghts_generate(small_lenghts, small_lenght_order);
            Wavefunction_relative_lenghts_generate(small_relative_lenghts, small_lenght_order); 
            }
        try
            {
            for (i = 0; i < count_orbitals; i++)
                {
                pointer_to_wavefunction = new T[wavefunction_size];
                pointer_to_probability = new T[wavefunction_size];
                pointer_to_Gradient = new T[wavefunction_size];
                pointers_to_wavefunctions[index[i]] = pointer_to_wavefunction;
                pointers_to_probabilities[index[i]] = pointer_to_probability;
                pointers_to_Gradients[index[i]] = pointer_to_Gradient;
                if (restriction == true and (spins[spin_paired[index[i]]] == -0.5)
                    and (bonding[index[i]] == -1 or Z[index[i]] == Z[bonding[index[i]]]))
                    {  // for restricted basis set method
                    pointers_to_wavefunctions[spin_paired[index[i]]] = pointer_to_wavefunction;
                    pointers_to_probabilities[spin_paired[index[i]]] = pointer_to_probability;
                    pointers_to_Gradients[spin_paired[index[i]]] = pointer_to_Gradient;
                    }
                }
            }
        catch (int)
            { // If allocation fail then deallocate allocated memory and return -1;
            for (i = 0; i < count_electrons; i++)
                {
                if (pointers_to_wavefunctions[i] != nullptr)
                    {
                    delete[] pointers_to_wavefunctions[i];
                    pointers_to_wavefunctions[i] = nullptr;
                    if (spin_paired[i] != -1 and bonding[i] == -1)
                        pointers_to_wavefunctions[spin_paired[i]] = nullptr;
                    }
                if (pointers_to_probabilities[i] != nullptr)
                    {
                    delete[] pointers_to_probabilities[i];
                    pointers_to_probabilities[i] = nullptr;
                    if (spin_paired[i] != -1 and bonding[i] == -1)
                        pointers_to_probabilities[spin_paired[i]] = nullptr;
                    }
                if (pointers_to_Gradients[i] != nullptr)
                    {
                    delete[] pointers_to_Gradients[i];
                    pointers_to_Gradients[i] = nullptr;
                    if (spin_paired[i] != -1 and bonding[i] == -1)
                        pointers_to_Gradients[spin_paired[i]] = nullptr;
                    }
                }
            return(-1);
            }
        try
            { // Allocating small wavefunctions for resonance integral lenght integration
            if (non_s1_system == true or s1_memory_optimization == false)
                for (i = 0; i < count_electrons; i++)
                    {
                    small_electron_numbers[i] = electron_numbers[i];
                    small_lenght_orders[i] = small_lenght_order;
                    small_x_range[i] = small_lenght_order;
                    small_y_range[i] = small_lenght_order;
                    small_z_range[i] = small_lenght_order;
                    small_wavefunctions[i] = new T[small_wavefunction_size];
                    small_probabilities[i] = new T[small_wavefunction_size];
                    small_atom_wavefunctions->electron_numbers.push_back(small_electron_numbers[i]);
                    small_atom_wavefunctions->lenght_orders.push_back(small_lenght_orders[i]);
                    small_atom_wavefunctions->x_range.push_back(small_x_range[i]);
                    small_atom_wavefunctions->y_range.push_back(small_y_range[i]);
                    small_atom_wavefunctions->z_range.push_back(small_z_range[i]);
                    small_atom_wavefunctions->wavefunctions.push_back(small_wavefunctions[i]);
                    small_atom_wavefunctions->probabilities.push_back(small_probabilities[i]);
                    small_atom_wavefunctions->x.push_back(atom_wavefunctions->x[i]);
                    small_atom_wavefunctions->y.push_back(atom_wavefunctions->y[i]);
                    small_atom_wavefunctions->z.push_back(atom_wavefunctions->z[i]);
                    small_atom_wavefunctions->n.push_back(atom_wavefunctions->n[i]);
                    small_atom_wavefunctions->l.push_back(atom_wavefunctions->l[i]);
                    small_atom_wavefunctions->m.push_back(atom_wavefunctions->m[i]);
                    small_atom_wavefunctions->Z.push_back(atom_wavefunctions->Z[i]);
                    }
                small_atom_wavefunctions_size = small_atom_wavefunctions->n.size();
            }
        catch (int)
            {
            for (i = 0; i < count_electrons; i++)
                {
                if (small_wavefunctions[i] != nullptr)
                    {
                    delete[] small_wavefunctions[i];
                    small_wavefunctions[i] = nullptr;
                    if (spin_paired[i] != -1 and bonding[i] == -1)
                        small_wavefunctions[spin_paired[i]] = nullptr;
                    }
                if (small_probabilities[i] != nullptr)
                    {
                    delete[] small_probabilities[i];
                    small_probabilities[i] = nullptr;
                    if (spin_paired[i] != -1 and bonding[i] == -1)
                        small_probabilities[spin_paired[i]] = nullptr;
                    }
                }
            return(-1);
            }
        }
    else
        {
        pointer_to_lenghts = atom_wavefunctions->lenghts[0];
        if (non_s1_system == true or s1_memory_optimization == false)
            small_lenghts = small_atom_wavefunctions->lenghts[0];
        for (i = 0; i < count_electrons; i++)
            {
            pointers_to_wavefunctions[i] = atom_wavefunctions->wavefunctions[i];
            pointers_to_probabilities[i] = atom_wavefunctions->probabilities[i];
            pointers_to_Gradients[i] = atom_wavefunctions->Gradients[i];
            }
        small_atom_wavefunctions_size = small_atom_wavefunctions->n.size();
        
        for (i = 0; i < small_atom_wavefunctions_size; i++)
            {
            small_electron_numbers[i] = small_atom_wavefunctions->electron_numbers[i];
            small_lenght_orders[i] = small_atom_wavefunctions->lenght_orders[i];
            small_x_range[i] = small_atom_wavefunctions->x_range[i];
            small_y_range[i] = small_atom_wavefunctions->y_range[i];
            small_z_range[i] = small_atom_wavefunctions->z_range[i];
            small_wavefunctions[i] = small_atom_wavefunctions->wavefunctions[i];
            small_probabilities[i] = small_atom_wavefunctions->probabilities[i];
            }
        }
    for (i = 0; i < count_orbitals; i++)
            {
            multipliers[i] = atom_wavefunctions->wavefunction_lenght_multipliers[index[i]];
            multiplier_array[i] = (vector_lenght/size_order) * multipliers[i];
            } // end of closed-shell basis set method optimalization code
    if (non_s1_system == true or s1_memory_optimization == false)
        for (i = 0; i < small_atom_wavefunctions_size; i++) // Generating small wavefunctions
            { 
            Orbitals_to_wavefunctions(small_atom_wavefunctions->n[i], small_atom_wavefunctions->l[i],
            small_atom_wavefunctions->m[i], small_lenght_orders[i], small_wavefunctions[i], small_lenghts,
            small_atom_wavefunctions->Z[i], atom_wavefunctions->wavefunction_lenght_multipliers[i] *
            (vector_lenght/small_lenght_order), &small_x_range[i], &small_y_range[i], &small_z_range[i]);
            Wavefunction_square(small_wavefunctions[i], small_probabilities[i], small_lenght_order);
            }
    // multithreading code
    #pragma omp parallel
        {
        #pragma omp for
        for (i = 0; i < count_orbitals; i++)
            {
            Orbitals_to_wavefunctions(n_i[i], l_i[i], m_i[i], size_order,
            pointers_to_wavefunctions[index[i]], pointer_to_lenghts, Z_i[i], multiplier_array[i],
            &x_range[index[i]], &y_range[index[i]], &z_range[index[i]]);
            }
        }
    if (compute_densities == true)
        {
        #pragma omp parallel
            {
            #pragma omp for
            for (i = 0; i < count_orbitals; i++)
                {
                Wavefunction_square(pointers_to_wavefunctions[index[i]], pointers_to_probabilities[index[i]],
                size_order);
                }
            }
        }
    if (allocate == true)
        {
        #pragma omp parallel
            {
            #pragma omp for
            for (i = 0; i < count_orbitals; i++)
                {
                Probabilities_thread(pointers_to_probabilities[index[i]],
                size_order, x_range[index[i]], y_range[index[i]], z_range[index[i]],
                &effective_radius_base_array[index[i]]);
                }
            }
        }
    #pragma omp parallel
        {
        #pragma omp for
        for (i = 0; i < count_orbitals; i++)
            {
            Gradient_thread(pointers_to_Gradients[index[i]], pointers_to_wavefunctions[index[i]], size_order);
            }
        }
    // end of multithreading code
    // closed-shell basis set method optimalization code    
    for (i = 0; i < count_orbitals; i++)
        if (restriction == true and (spins[spin_paired[index[i]]] == -0.5)
        and (bonding[index[i]] == -1 or Z[index[i]] == Z[bonding[index[i]]])) // for restricted basis set method
            {
            effective_radius_base_array[spin_paired[index[i]]] = effective_radius_base_array[index[i]];
            x_range[spin_paired[index[i]]] = x_range[index[i]];
            y_range[spin_paired[index[i]]] = y_range[index[i]];
            z_range[spin_paired[index[i]]] = z_range[index[i]];
            if (non_s1_system == true or s1_memory_optimization == false)
                {
                small_x_range[spin_paired[index[i]]] = small_x_range[index[i]];
                small_y_range[spin_paired[index[i]]] = small_y_range[index[i]];
                small_z_range[spin_paired[index[i]]] = small_z_range[index[i]];
                }
            }
    // end of closed-shell basis set method optimalization code
    if (allocate == true)
        {
        for (i = 0; i < count_electrons; i++) // copy values to vectors for restricted and unrestricted method
            {
            atom_wavefunctions->wavefunctions.push_back(pointers_to_wavefunctions[i]);
            atom_wavefunctions->probabilities.push_back(pointers_to_probabilities[i]);
            atom_wavefunctions->effective_radius_base.push_back(effective_radius_base_array[i]);
            atom_wavefunctions->Gradients.push_back(pointers_to_Gradients[i]);
            }
        }
    return(0);
    }
// end of section 5 - generating matrices of integrals and Fock matrices
//section 6: completing basis set method and user interface handling
template <typename T>
int Slater_basis_set_calculations<T>::Atom_orbitals_generate(string UI_input, atom_orbitals *atom_orbitals_PTR)
    { // convert strings of atoms, ions and radicals
    unsigned int i;
    unsigned int input_size;
    unsigned int element_number;
    
    int shift, sign, string_position;
    vector<int> positions_to_erase;
    string input, shift_string;
    bool found;
    
    struct numbers {
    vector<unsigned int> numbers;
    string text;
    } numbers;
    numbers.numbers = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    numbers.text = {"0123456789"};
    
    struct signs {
    vector<int> sign;
    string text;
    } signs;
    signs.sign = { 1, -1};
    signs.text = {"+-"};
    
    struct order_elements {
    vector<unsigned int> numbers;
    vector<string> elements;
    } order_elements;
    
    auto reserve = [&atom_orbitals_PTR](unsigned int size){
    atom_orbitals_PTR->n.reserve(size);
    atom_orbitals_PTR->l.reserve(size);
    atom_orbitals_PTR->m.reserve(size);
    atom_orbitals_PTR->s.reserve(size);
    atom_orbitals_PTR->wavefunction_lenght_multipliers.reserve(size);
    atom_orbitals_PTR->wavefunction_coefficients.reserve(size);
    atom_orbitals_PTR->bonding.reserve(size);
    atom_orbitals_PTR->paired_with_previous.reserve(size);
    return(0);
    };
    
    order_elements.numbers = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
    81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108,
    109, 110, 111, 112, 113, 114, 115, 116, 117};
    order_elements.elements = {"H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
    
    found = false;
    shift = 1;
    shift_string = "";
    sign = 0;
    element_number = 0;
    input_size = UI_input.size();
    input = UI_input;
    
    for (i = 0; i < input_size; i++) // search for ionts signatures
        {
        string_position = numbers.text.find(input[i]); // serch for numbers
        if (string_position != -1)
            {
            shift_string = shift_string + numbers.text[string_position];
            positions_to_erase.push_back(i);
            }
        string_position = signs.text.find(input[i]); // serch for signs
        if (string_position != -1)
            {
            sign = signs.sign[string_position];
            positions_to_erase.push_back(i);
            }
        }
    if (shift_string.size() > 0)
        shift = stoi(shift_string);
    else
        shift = 1;
    
    for (i = positions_to_erase.size(); i > 0; i--) // removing numbers and signs from input string
        input.erase(positions_to_erase[i - 1], 1);
        
    for (i = 0; i < order_elements.elements.size(); i++) // finding electron configuration of atoms and ions
        if ((found == false) and (order_elements.elements[i] == input))
            {
            element_number = order_elements.numbers[i];
            found = true;
            }
    if (found == false)
        return(-1);
        
    shift = shift * -sign; // detection of wrong shifts
    if (element_number + shift < 0 or element_number + shift >= 118) // not electrons or not discovered element
        return(-1);
        
    if ((element_number < 18 and (element_number + shift) > 18) or // exceeding of rare gas configuration
    (element_number < 18 and (element_number + shift) > 18) or
    (element_number < 36 and (element_number + shift) > 36) or
    (element_number < 54 and (element_number + shift) > 54) or
    (element_number < 86 and (element_number + shift) > 86))
        return(-1);
    
    if (element_number > 10 and vector_lenght < 10) // Simulations for large space
        vector_lenght = 10;
    if (element_number > 36 and vector_lenght < 15)
        vector_lenght = 15;
    
    input.clear();
    input = order_elements.elements[element_number + shift];
    reserve(element_number);
    Atoms_to_valence_orbitals(input, atom_orbitals_PTR);
    atom_orbitals_PTR->charge = -shift;
    atom_orbitals_PTR->Z = element_number + 1;
    atom_orbitals_PTR->reduced_Z = atom_orbitals_PTR->reduced_Z - shift;
    return(0);
    }
template <typename T>
T Slater_basis_set_calculations<T>::Nucleus_repulsive_energy(atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j;
    unsigned int sum_electrons;
    vector<unsigned int> ind_nuc;
    T square_lenght;
    T energy;
    T* x = atom_wavefunctions->x.data();
    T* y = atom_wavefunctions->y.data();
    T* z = atom_wavefunctions->z.data();
    unsigned int* Z = atom_wavefunctions->reduced_Z.data();
    
    sum_electrons = atom_wavefunctions->n.size();
    energy = 0;
    i = 0;
    ind_nuc.push_back(0);
    
    for (i = 1; i < sum_electrons; i++)
        if (x[i] != x[i-1] or y[i] != y[i-1] or z[i] != z[i-1])
            ind_nuc.push_back(i);
      
    for (i = 0; i < ind_nuc.size(); i++)
        for (j = i + 1; j < ind_nuc.size(); j++)
            {
            if (j != i and ((x[ind_nuc[j]] != x[ind_nuc[i]]) or (y[ind_nuc[j]] != y[ind_nuc[i]])
            or (z[ind_nuc[j]] - z[ind_nuc[i]])))
                {
                square_lenght = (x[ind_nuc[j]] - x[ind_nuc[i]]) * (x[ind_nuc[j]] - x[ind_nuc[i]])
                + (y[ind_nuc[j]] - y[ind_nuc[i]]) * (y[ind_nuc[j]] - y[ind_nuc[i]])
                + (z[ind_nuc[j]] - z[ind_nuc[i]]) * (z[ind_nuc[j]] - z[ind_nuc[i]]);
                square_lenght = sqrt(abs(square_lenght)) * Hartree_lenght;
                energy = energy + (e*e/(4*Pi*E0) * Z[ind_nuc[i]] * (Z[ind_nuc[j]]/square_lenght));
                }
            }
    return(energy);
    }
template <typename T>
int Slater_basis_set_calculations<T>::String_to_list_electrons(string UI_input, unsigned int size_order,
bool extern_coordinates, vector<T>* x_2, vector<T>* y_2, vector<T>* z_2)
    {
    unsigned int i, j, k;
    unsigned int last_readed_index;
    unsigned int input_size;
    unsigned int matrix_order;
    unsigned int count_atoms, count_coordinates, count_bonds, count_potentials;
    
    bool read_switch;
    bool previous_deleted;
    
    string input;
    string character;
    array<string, max_atoms> atoms_string;
    array<string, max_atoms> coordinates_string;
    array<string, max_atoms> bonds_string;
    array<string, max_atoms> potentials_string;
    
    vector<atom_orbitals> chain;
    vector<atom_wavefunctions> wavefunctions;
    vector<T> x, y, z;
    vector<T> potentials;
    
    const string characters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789+-";
    const string numbers = "+-0123456789.";
    const string begin_brackets = "([{<";
    const string end_brackets = ")]}>";
    
    struct bond { // for molecules
    vector<unsigned int> atom_1;
    vector<unsigned int> atom_2;
    vector<unsigned int> bond_count;
    } bond_list;
    
    count_atoms = 0;
    count_coordinates = 0;
    count_bonds = 0;
    count_potentials = 0;
    read_switch = true;
    character = " ";
    matrix_order = 0;
    input_size = UI_input.size();
    input = UI_input;
        
    if (input_size <= 3 or ((input_size <= 4) and (input[input_size -1] == '+' or input[input_size -1] == '-')))// For one atom or ion
        if (Atom_orbitals_generate(input, &atoms) != -1)
            {
            Create_atomic_wavefunctions(&atoms, &results, size_order, 0, 0, 0);
            matrix_order = atoms.n.size();
            }
        else
            return(-1);
    else // for molecules
        {
        last_readed_index = 10000;
        for (i = 0; i < input_size; i++) // reading list of atoms
            {
            for (j = 0; j < characters.size(); j++) // check for UI_input[i] in list of characters
                if (UI_input[i] == characters[j])
                    if (read_switch == true)
                        {
                        if ((i - last_readed_index) != 1)
                            count_atoms++;
                            
                        character[0] = UI_input[i];
                        atoms_string[count_atoms - 1].append(character);
                        last_readed_index = i;
                        }
            for (j = 0; j < begin_brackets.size(); j++) // check for UI_input[i] in begin_brackets
                if (UI_input[i] == begin_brackets[j])
                    read_switch = false;
            for (j = 0; j < begin_brackets.size(); j++) // check for UI_input[i] in end_brackets
                if (UI_input[i] == end_brackets[j])
                    read_switch = true;
            }
        if (extern_coordinates == false)
            {
            read_switch = false;
            last_readed_index = 10000;
            for (i = 0; i < input_size; i++) // reading list of coordinates
                {
                for (j = 0; j < numbers.size(); j++) // check for UI_input[i] in list of numbers
                    if (UI_input[i] == numbers[j] and read_switch == true)
                        {
                        if ((i - last_readed_index) != 1)
                            count_coordinates++;
                            
                        character[0] = UI_input[i];
                        coordinates_string[count_coordinates - 1].append(character);
                        last_readed_index = i;
                        }
                if (UI_input[i] == begin_brackets[1])
                    read_switch = true;
                if (UI_input[i] == end_brackets[1])
                    read_switch = false;
                }
            }
        read_switch = false;
        last_readed_index = 10000;
        for (i = 0; i < input_size; i++) // reading list of bounds
            {
            for (j = 0; j < numbers.size(); j++) // check for UI_input[i] in list of numbers
                if (UI_input[i] == numbers[j] and read_switch == true)
                    {
                    if ((i - last_readed_index) != 1)
                        count_bonds++;
                        
                    character[0] = UI_input[i];
                    bonds_string[count_bonds - 1].append(character);
                    last_readed_index = i;
                    }
            if (UI_input[i] == begin_brackets[2])
                read_switch = true;
            if (UI_input[i] == end_brackets[2])
                read_switch = false;
            }
        read_switch = false;
        last_readed_index = 10000;
        for (i = 0; i < input_size; i++) // reading list of potentials
            {
            for (j = 0; j < numbers.size(); j++) // check for UI_input[i] in list of numbers
                if (UI_input[i] == numbers[j] and read_switch == true)
                    {
                    if ((i - last_readed_index) != 1)
                        count_potentials++;
                        
                    character[0] = UI_input[i];
                    potentials_string[count_potentials - 1].append(character);
                    last_readed_index = i;
                    }
            if (UI_input[i] == begin_brackets[0])
                read_switch = true;
            if (UI_input[i] == end_brackets[0])
                read_switch = false;
            }
        chain.reserve(count_atoms); // reserve necessary amount of memory to vectors;
        wavefunctions.resize(count_atoms);
        
        if (extern_coordinates == true) // copy coordinates from extern source
            for (i = 0; i < x_2->size(); i++)
                {
                x.push_back(x_2->operator[](i));
                y.push_back(y_2->operator[](i));
                z.push_back(z_2->operator[](i));
                }
        else
            {
            x.push_back(0);
            y.push_back(0);
            z.push_back(0);
            for (i = 0; i < count_coordinates; i++) // convert coordinates to double
                {
                if (i % 3 == 0)
                    x.push_back(stod(coordinates_string[i]));
                if (i % 3 == 1)
                    y.push_back(stod(coordinates_string[i]));
                if (i % 3 == 2)
                    z.push_back(stod(coordinates_string[i]));
                }
                if (count_coordinates % 3 == 1)
                    {
                    y.push_back(0);
                    z.push_back(0);
                    }
                if (count_coordinates % 3 == 2)
                    z.push_back(0);
            }  
        for (i = 0; i < count_atoms; i++) // convert list of atoms_string to list of wavefunctions and applying coordinates
            if (Atom_orbitals_generate(atoms_string[i], &chain[i]) != -1)
                {
                Create_atomic_wavefunctions(&chain[i], &wavefunctions[i], size_order, x[i], y[i], z[i]);
                matrix_order = matrix_order + wavefunctions[i].n.size();
                }
            else
                if (atoms_string[i] == "H+")
                    {
                    Atom_orbitals_generate("H", &chain[i]);
                    Create_atomic_wavefunctions(&chain[i], &wavefunctions[i], size_order, x[i], y[i], z[i]);
                    wavefunctions[i].wavefunction_coefficients[0] = 0.00;
                    wavefunctions[i].wavefunction_constraints[0] = 1;
                    matrix_order = matrix_order + wavefunctions[i].n.size();
                    }
                else
                    return(-1); // check for right input atoms strings
        if (count_bonds >= 1)
            bonded_system = true;
        for (i = 0; i < count_bonds; i++) // convert bonds to list of bonds
            {
            if (i % 3 == 0)
                bond_list.atom_1.push_back(stoi(bonds_string[i]));
            if (i % 3 == 1)
                bond_list.atom_2.push_back(stoi(bonds_string[i]));
            if (i % 3 == 2)
                bond_list.bond_count.push_back(stoi(bonds_string[i]));
            }
        for (i = 0; i < count_potentials; i++) // convert potentials to double
            potentials.push_back(stod(potentials_string[i]));
        
        if ((x.size()  < chain.size()) or (y.size()  < chain.size()) or (z.size() < chain.size()))  
            return(-1); // check for a sufficient number of coordinates
            
        for (i = 0; i < bond_list.atom_1.size(); i++) // applying bonds
            for (j = 0; j < bond_list.bond_count[i]; j++)
                {
                Create_bond_atomic_wavefunctions(&wavefunctions[bond_list.atom_1[i] - 1], &wavefunctions[bond_list.atom_2[i] - 1],
                1, chain[bond_list.atom_1[i] - 1].electronegativity, chain[bond_list.atom_2[i] - 1].electronegativity,
                wavefunctions[bond_list.atom_1[i] - 1].x[0] - wavefunctions[bond_list.atom_2[i] - 1].x[0],
                wavefunctions[bond_list.atom_1[i] - 1].y[0] - wavefunctions[bond_list.atom_2[i] - 1].y[0],
                wavefunctions[bond_list.atom_1[i] - 1].z[0] - wavefunctions[bond_list.atom_2[i] - 1].z[0]);
                }
        // reserve memory for results
        results.lenghts.reserve(1);
        results.wavefunctions.reserve(electron_number);
        results.probabilities.reserve(electron_number);
        results.Gradients.reserve(electron_number);
        results.lenght_orders.reserve(electron_number);
        results.x_range.reserve(electron_number);
        results.y_range.reserve(electron_number);
        results.z_range.reserve(electron_number);
        results.wavefunction_coefficients.reserve(electron_number);
        results.wavefunction_lenght_multipliers.reserve(electron_number);
        results.effective_radius_base.reserve(electron_number);
        results.spins.reserve(electron_number);
        results.spin_paired.reserve(electron_number);
        results.bonding.reserve(electron_number);
        results.antibonding.reserve(electron_number);
        results.pi_bonding.reserve(electron_number);
        results.wavefunction_constraints.reserve(electron_number);
        results.n.reserve(electron_number);
        results.l.reserve(electron_number);
        results.m.reserve(electron_number);
        results.charge.reserve(electron_number);
        results.reduced_Z.reserve(electron_number);
        results.Z.reserve(electron_number);
        results.electron_numbers.reserve(electron_number);
        results.x.reserve(electron_number);
        results.y.reserve(electron_number);
        results.z.reserve(electron_number);
        for (i = 0; i < count_atoms; i++) // sumarizing wavefunctions
            Sum_atomic_wavefunctions(&results, &wavefunctions[i]);
        }
    try {
        nuclear_atraction_integral_matrix = new T[matrix_order * matrix_order];
        nucleuses_atractions = new T[matrix_order * matrix_order];
        nucleuses_distances = new T[matrix_order * matrix_order];
        coulombic_integral_matrix = new T[matrix_order * matrix_order];
        overlap_integral_matrix = new T[matrix_order * matrix_order];
        overlap_effective_lenght_integral_matrix = new T[matrix_order * matrix_order];
        resonance_integral_matrix = new T[matrix_order * matrix_order];
        kinetic_integral_matrix = new T[matrix_order * matrix_order];
        basis_set_matrix = new T[matrix_order * matrix_order];
        correction_matrix = new T[matrix_order * matrix_order];
        corr_basis_set_matrix = new T[matrix_order * matrix_order];
        spin_density_matrix = new T[matrix_order * matrix_order];
        memset(resonance_integral_matrix, 0, matrix_order * matrix_order);
        }
    catch (int)
        {
        if (nuclear_atraction_integral_matrix != nullptr){
            delete[] nuclear_atraction_integral_matrix;
            nuclear_atraction_integral_matrix = nullptr;
            }
        if (nucleuses_distances != nullptr){
            delete[] nucleuses_distances;
            nucleuses_distances = nullptr;
            }
        if (nucleuses_atractions != nullptr){
            delete[] nucleuses_atractions;
            nucleuses_atractions = nullptr;
            }
        if (coulombic_integral_matrix != nullptr){
            delete[] coulombic_integral_matrix;
            coulombic_integral_matrix = nullptr;
            }
        if (overlap_integral_matrix != nullptr){
            delete[] overlap_integral_matrix;
            overlap_integral_matrix = nullptr;
            }
        if (overlap_effective_lenght_integral_matrix != nullptr){
            delete[] overlap_effective_lenght_integral_matrix;
            overlap_effective_lenght_integral_matrix = nullptr;
            }
        if (resonance_integral_matrix != nullptr){
            delete[] resonance_integral_matrix;
            resonance_integral_matrix = nullptr;
            }
        if (kinetic_integral_matrix != nullptr){
            delete[] kinetic_integral_matrix;
            kinetic_integral_matrix = nullptr;
            }
        if (basis_set_matrix != nullptr){
            delete[] basis_set_matrix;
            basis_set_matrix = nullptr;
            }
        if (correction_matrix != nullptr){
            delete[] correction_matrix;
            correction_matrix = nullptr;
            }
        if (corr_basis_set_matrix != nullptr){
            delete[] corr_basis_set_matrix;
            corr_basis_set_matrix = nullptr;
            }
        if (spin_density_matrix != nullptr){
            delete[] spin_density_matrix;
            spin_density_matrix = nullptr;
            }
        return(-1);
        }
    j = 0;
    k = 0;
    for (i = 0; i < (matrix_order * matrix_order); i++)
        correction_matrix[i] = 0;
    
    for (i = 0; i < potentials.size(); i++) // creating corr_basis_set_matrix
        {
        while (k < wavefunctions[i].count_electrons[0])
            {
            correction_matrix[(j + k) * (matrix_order + 1)] = potentials[i] * e;
            k++;
            }
        j = j + k;
        k = 0;
        }
    return(0);
    }
template <typename T>
T Slater_basis_set_calculations<T>::Calculate(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool allocate, vector<T>* values)
    {
    unsigned int i, j;
    unsigned int matrix_order;
    T Hamiltonian, last_Hamiltonian;
    
    last_Hamiltonian = 0;
    matrix_order = results.n.size();
    memset(basis_set_matrix, 0, matrix_order * matrix_order);
    memset(corr_basis_set_matrix, 0, matrix_order * matrix_order);
    memset(spin_density_matrix, 0, matrix_order * matrix_order);
    if (allocate == true)
        if (Generate_atomic_wavefunctions(&results, &small_results, size_order, true, true) == -1)
            return(-1);  
        else
            {
            Calculate_basis_set_matrix(nuclear_atraction_integral_matrix, coulombic_integral_matrix,
            resonance_integral_matrix,
            kinetic_integral_matrix, basis_set_matrix, matrix_order, &results);
            Calculate_corr_basis_set_matrix(basis_set_matrix, correction_matrix, corr_basis_set_matrix, matrix_order);
            Hamiltonian = Solve_basis_set_matrix(corr_basis_set_matrix, overlap_integral_matrix, matrix_order,
            values, &results);
            nucleus_repulsive_energy = Nucleus_repulsive_energy(&results);
            Hamiltonian = Hamiltonian + nucleus_repulsive_energy;
            }
    for (i = 0; i < max_iterations + 1; i++)
        {
        iterations++;
        Create_nuclear_atraction_integral_matrix(nuclear_atraction_integral_matrix, nucleuses_atractions,
        matrix_order, &results);
        Create_coulombic_integral_matrix(coulombic_integral_matrix, matrix_order, &results, &small_results);
        if (bonded_system == true)
            {
            Create_overlap_integral_matrix(overlap_integral_matrix, matrix_order, &results);
            Calculate_resonance_integral_matrix(overlap_integral_matrix, overlap_effective_lenght_integral_matrix,
            resonance_integral_matrix, matrix_order, &results, &small_results);
            }
        Calculate_kinetic_integral_matrix(kinetic_integral_matrix, matrix_order, &results);
        
        Calculate_basis_set_matrix(nuclear_atraction_integral_matrix, coulombic_integral_matrix,
        resonance_integral_matrix, kinetic_integral_matrix, basis_set_matrix, matrix_order, &results);
        Calculate_corr_basis_set_matrix(basis_set_matrix, correction_matrix, corr_basis_set_matrix, matrix_order);
        Hamiltonian = Solve_basis_set_matrix(corr_basis_set_matrix, overlap_integral_matrix,
        matrix_order, values, &results);
        if ((abs(Hamiltonian/last_Hamiltonian) < (2 - minimal_fidelity)) and
        (abs(Hamiltonian/last_Hamiltonian) > minimal_fidelity)
        and (Hamiltonian/last_Hamiltonian > 0))
            break;
        
        last_Hamiltonian = Hamiltonian;
        Generate_atomic_wavefunctions(&results, &small_results, size_order, false, true);
        }
    nucleus_repulsive_energy = Nucleus_repulsive_energy(&results);
    Hamiltonian = Hamiltonian + nucleus_repulsive_energy;
    return(Hamiltonian);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Create_excitation(unsigned int electron_number,
unsigned int n, unsigned int l, int m, T spin, bool generate)
    {
    unsigned int i;
    unsigned int x_range, y_range, z_range;
    unsigned int size = results.n.size();
    bool restriction = true;
    bool non_s1 = false;
    T effective_radius_base;
    
    // Check for spin-pairing
    for (i = 0; i < size; i++)
        {
        if (results.spin_paired[i] == -1)
            restriction = false;
        if (results.n[i] > 1)
            non_s1 = true;
        }
    
    // searching electron
    for (i = 0; i < size; i++)
        if (results.electron_numbers[i] == electron_number)
            break;
    
    // avoid generation excitation for closed-shell grid with allocated memory
    if (restriction == true)
        {
        if (iterations > 0)
            return(-1);
        else
           results.spin_paired[i] = -1;
        }
    // avoid generation excitation from s1 electron systems with not allocated for small_results grid
    if (non_s1 == false)
        {
        // before iterations set a class switch s1_memory-Optimization to false
        if (iterations > 0 and s1_memory_optimization == true)
            return(-1);
        else
            s1_memory_optimization = false;
        }
    // Set new set of quantum numbers
    results.n[i] = n;
    results.l[i] = l;
    results.m[i] = m;
    results.spins[i] = spin;
    
    results.wavefunction_lenght_multipliers[i] = 1;
    if (generate == true)
        {
        // Generate a new wavefunction
        results.wavefunction_lenght_multipliers[i] = 1;
        Orbitals_to_wavefunctions(n, l, m, results.lenght_orders[i],
        results.wavefunctions[i], results.lenghts[0], results.Z[i], 1, &x_range, &y_range, &z_range);
        results.x_range[i] = x_range;
        results.y_range[i] = y_range;
        results.z_range[i] = z_range;
        Wavefunction_range_detect(small_results.wavefunctions[i], small_results.lenght_orders[i],
        &x_range, &y_range, &z_range);
        small_results.x_range[i] = x_range;
        small_results.y_range[i] = y_range;
        small_results.z_range[i] = z_range;
        results.x_range[i] = x_range;
        results.y_range[i] = y_range;
        results.z_range[i] = z_range;
        // Generate new electron density
        Wavefunction_square(results.wavefunctions[i], results.probabilities[i], results.lenght_orders[i]);
        // Recompute a effective radius_base
        Probabilities_thread(results.probabilities[i], results.lenght_orders[i],
        results.x_range[i], results.y_range[i], results.z_range[i], &effective_radius_base);
        results.effective_radius_base[i] = effective_radius_base;
        // Recompute a gradient
        Gradient_thread(results.Gradients[i], results.wavefunctions[i], results.lenght_orders[i]);
        
        }
    return(0);
    }
template <typename T>
int Slater_basis_set_calculations<T>::Clear()
    {
    unsigned int i, j;
    unsigned int order = results.n.size();
    bool previous_deleted;
    
    for (i = 0; i < results.wavefunctions.size(); i++) {
            previous_deleted = false; // avoid double delete
            for (j = 0; j < i; j++)
                if (results.wavefunctions[i] == results.wavefunctions[j])
                    previous_deleted = true;
            
            if (previous_deleted == false) {
                delete[] results.wavefunctions[i];
                delete[] results.probabilities[i];
                delete[] results.Gradients[i];
                }
            }
    for (i = 0; i < results.lenghts.size(); i++)
        delete[] results.lenghts[i];
    
    for (i = 0; i < small_results.lenghts.size(); i++) {
        if (small_results.lenghts[i] != nullptr)
            delete[] small_results.lenghts[i];
        }
    for (i = 0; i < small_results.relative_lenghts.size(); i++) {
        if (small_results.relative_lenghts[i] != nullptr)
            delete[] small_results.relative_lenghts[i];
        }
    for (i = 0; i < small_results.wavefunctions.size(); i++) {
        if (small_results.wavefunctions[i] != nullptr)
            delete[] small_results.wavefunctions[i];
        if (small_results.probabilities[i] != nullptr)
            delete[] small_results.probabilities[i];
        }
    if (nuclear_atraction_integral_matrix != nullptr) {
        delete[] nuclear_atraction_integral_matrix;
        nuclear_atraction_integral_matrix = nullptr;
        }
    if (nucleuses_distances != nullptr) {
        delete[] nucleuses_distances;
        nucleuses_distances = nullptr;
        }
    if (nucleuses_atractions != nullptr) {
        delete[] nucleuses_atractions;
        nucleuses_atractions = nullptr;
        }
    if (coulombic_integral_matrix != nullptr) {
        delete[] coulombic_integral_matrix;
        coulombic_integral_matrix = nullptr;
        }
    if (overlap_integral_matrix != nullptr) {
        delete[] overlap_integral_matrix;
        overlap_integral_matrix = nullptr;
        }
    if (overlap_effective_lenght_integral_matrix != nullptr) {
        delete[] overlap_effective_lenght_integral_matrix;
        overlap_effective_lenght_integral_matrix = nullptr;
        }
    if (resonance_integral_matrix != nullptr) {
        delete[] resonance_integral_matrix;
        resonance_integral_matrix = nullptr;
        }
    if (kinetic_integral_matrix != nullptr) {
        delete[] kinetic_integral_matrix;
        kinetic_integral_matrix = nullptr;
        }
    if (basis_set_matrix != nullptr) {
        delete[] basis_set_matrix;
        basis_set_matrix = nullptr;
        }
    if (correction_matrix != nullptr) {
        delete[] correction_matrix;
        correction_matrix = nullptr;
        }
    if (corr_basis_set_matrix != nullptr) {
        delete[] corr_basis_set_matrix;
        corr_basis_set_matrix = nullptr;
        }
    if (spin_density_matrix != nullptr) {
        delete[] spin_density_matrix;
        spin_density_matrix = nullptr;
        }
    determinants.clear();
    spectra_EPR.clear();
    electron_spectra.clear();
    
    atoms.n.clear();
    atoms.l.clear();
    atoms.m.clear();
    atoms.s.clear();
    atoms.wavefunction_lenght_multipliers.clear();
    atoms.wavefunction_coefficients.clear();
    atoms.bonding.clear();
    atoms.paired_with_previous.clear();
    
    results.lenghts.clear();
    results.wavefunctions.clear();
    results.probabilities.clear();
    results.Gradients.clear();
    results.lenght_orders.clear();
    results.x_range.clear();
    results.y_range.clear();
    results.z_range.clear();
    results.wavefunction_coefficients.clear();
    results.wavefunction_lenght_multipliers.clear();
    results.effective_radius_base.clear();
    results.spins.clear();
    results.spin_paired.clear();
    results.bonding.clear();
    results.antibonding.clear();
    results.pi_bonding.clear();
    results.wavefunction_constraints.clear();
    results.charge.clear();
    results.n.clear();
    results.l.clear();
    results.m.clear();
    results.count_electrons.clear();
    results.reduced_Z.clear();
    results.Z.clear();
    results.electron_numbers.clear();
    results.x.clear();
    results.y.clear();
    results.z.clear();
    
    
    small_results.electron_numbers.clear();
    small_results.lenght_orders.clear();
    small_results.x_range.clear();
    small_results.y_range.clear();
    small_results.z_range.clear();
    small_results.lenghts.clear();
    small_results.relative_lenghts.clear();
    small_results.wavefunctions.clear();
    small_results.probabilities.clear();
    small_results.x.clear();
    small_results.y.clear();
    small_results.z.clear();
    small_results.n.clear();
    small_results.l.clear();
    small_results.m.clear();
    small_results.Z.clear();
    electron_number = 0;
    iterations = 0;
    determinant_exception_handle = 0;
    nucleus_repulsive_energy = 0;
    relative_permitivity = 1;
    return(0);
    }
template <typename T>
Slater_basis_set_calculations<T>::~Slater_basis_set_calculations(){
Clear();}
template class Slater_basis_set_calculations<double>; /*
Author of this source code Ing. Pavel Florian Ph.D. licensed this source code under the the Apache License:
Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/ */
