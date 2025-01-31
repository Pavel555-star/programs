#include "basis_set_calculations.h"
using namespace std;
template <typename T>
inline T basis_set_calculations<T>::Determinant(unsigned int order, T* pointer, T* buffer, T* denominator, T* temp1, T* temp2) 
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
T basis_set_calculations<T>::basis_set_Determinant_set(unsigned int order, T* pointer,unsigned int count, T min, T step,
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
T basis_set_calculations<T>::basis_set_Determinant_solver(unsigned int order, T* pointer)
    {   // multithreading code
        unsigned int i;
        T* determinant_values = nullptr;
        
        try {
            determinant_values = new T[10000 * order];
            } // Calculating of determinant values array
        catch (int){
            determinant_exception_handle = 1;
            return (-1);
            }
        thread t136(&basis_set_calculations::basis_set_Determinant_set, this, order, pointer, 1250 * order, -1.000 * order,
        0.0002, determinant_values);
        thread t137(&basis_set_calculations::basis_set_Determinant_set, this, order, pointer, 1250 * order, -0.750 * order,
        0.0002, determinant_values + 1250 * order);
        thread t138(&basis_set_calculations::basis_set_Determinant_set, this, order, pointer, 1250 * order, -0.500 * order,
        0.0002, determinant_values + 2500 * order);
        thread t139(&basis_set_calculations::basis_set_Determinant_set, this, order, pointer, 1250 * order, -0.250 * order,
        0.0002, determinant_values + 3750 * order);
        thread t140(&basis_set_calculations::basis_set_Determinant_set, this, order, pointer, 1250 * order,  0.000 * order,
        0.0002, determinant_values + 5000 * order);
        thread t141(&basis_set_calculations::basis_set_Determinant_set, this, order, pointer, 1250 * order,  0.250 * order,
        0.0002, determinant_values + 6250 * order);
        thread t142(&basis_set_calculations::basis_set_Determinant_set, this, order, pointer, 1250 * order,  0.500 * order,
        0.0002, determinant_values + 7500 * order);
        thread t143(&basis_set_calculations::basis_set_Determinant_set, this, order, pointer, 1250 * order,  0.750 * order,
        0.0002, determinant_values + 8750 * order);
        t136.join();
        t137.join();
        t138.join();
        t139.join();
        t140.join();
        t141.join();
        t142.join();
        t143.join();
        // end of multithreading code
        determinants.clear();
        for (i =1; i < 10000 * order; i++)
            {
            if ((determinant_values[i] * determinant_values[i] <= 0.01) and (determinant_values[i] * determinant_values[i] <
            determinant_values[i + 1] * determinant_values[i + 1]) and (determinant_values[i] * determinant_values[i] <
            determinant_values[i - 1] * determinant_values[i - 1]))
                determinants.push_back((double(i)/5000) -1 * order);
            }
        delete[] determinant_values;
        return(0);
    }
// End of Section 1 - work fith Fock matrix - inherited from Huckel_calculations
// Section 2 - generating the wavefunctions
template <typename T>
T basis_set_calculations<T>::Wavefunction_lenghts_generate(T* lenghts, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_1s_generate(T* wavefunction, T* lenghts, int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_2s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_3s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7s_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_2px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_2pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_2py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_3px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_3pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_3py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7px_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7pz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7py_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_3dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_3dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_3dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_3dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_3dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7dx2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7dxz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7dz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7dyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7dxy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_4fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_6fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7fx_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7fz_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7fxz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7fz3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7fyz2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7fxyz_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_7fy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier,
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
T basis_set_calculations<T>::Wavefunction_5gz4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
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
T basis_set_calculations<T>::Wavefunction_5gz3y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
                x_i = (int(x) - int(lenght_order));
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
T basis_set_calculations<T>::Wavefunction_5gz3x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
            y_i = (int(y) - int(lenght_order));
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
T basis_set_calculations<T>::Wavefunction_5gz2xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_5gz2_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_5gz_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_5gz_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_5gxy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_5gx4_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_6gz4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
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
T basis_set_calculations<T>::Wavefunction_6gz3y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
                x_i = (int(x) - int(lenght_order));
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
T basis_set_calculations<T>::Wavefunction_6gz3x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
            y_i = (int(y) - int(lenght_order));
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
T basis_set_calculations<T>::Wavefunction_6gz2xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_6gz2_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_6gz_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_6gz_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_6gxy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_6gx4_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_7gz4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
            y_i = (int(y) - int(lenght_order));
            for (x = 0; x < side; x++)
                {
                x_i = (int(x) - int(lenght_order));
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
T basis_set_calculations<T>::Wavefunction_7gz3y_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
                x_i = (int(x) - int(lenght_order));
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
T basis_set_calculations<T>::Wavefunction_7gz3x_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
            y_i = (int(y) - int(lenght_order));
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
T basis_set_calculations<T>::Wavefunction_7gz2xy_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_7gz2_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_7gz_y3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_7gz_x3_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_7gxy_x2_y2_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_7gx4_y4_generate(T* wavefunction, T* lenghts, unsigned int Z, T multiplier, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_normalize(T* wavefunction_pointer, T normalisation_constant, unsigned int size)
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
T basis_set_calculations<T>::Orbitals_to_wavefunctions(unsigned int n, unsigned int l, int m, unsigned int lenght_order,
T* wavefunction, T* lenghts, unsigned int Z, T multiplier)
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
        }
    size = (2 * lenght_order + 1) * (2 * lenght_order + 1) * (2 * lenght_order + 1);
    Wavefunction_normalize(wavefunction, normalisation_constant, size);
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Get_relative_Hartree_length(unsigned int Z, unsigned int n)
    { // Setting relative permeability and permeability constant for wavefunction calculations
    T relative_Hartree_lenght;
    T relative_electron_mass;
    T relativistic_shrinkage;
    
    relative_electron_mass = (Z * mp)/(Z * mp + me);
    relativistic_shrinkage = 1/sqrt(1 - pow(Z * hyperfine_structure_constant/n, 2));
    relative_Hartree_lenght = relativistic_shrinkage/relative_electron_mass;
    return(relative_Hartree_lenght);
    }
// End of Section 2 - generating wavefunctions, Section 3 - mathematical operations for wavefunctions, probabilities densities and
// integrals
template <typename T>
unsigned int basis_set_calculations<T>::Wavefunction_multiply(T* wavefunction_1, T* wavefunction_2, T* probabilities,
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
T basis_set_calculations<T>::Wavefunction_multiply(T* wavefunction_1, T* wavefunction_2, T* probabilities, unsigned int lenght_order)
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
T basis_set_calculations<T>::Wavefunction_relative_lenghts_generate(T* reverse_relative_lenghts,
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
                    pre_index = (i * lenght_order * lenght_order * lenght_order * lenght_order * lenght_order) +
                    (j * lenght_order * lenght_order * lenght_order * lenght_order) + (k * lenght_order * lenght_order * lenght_order) +
                    (l * lenght_order * lenght_order * lenght_order * lenght_order);
                    for (m = 0; m < side; m++)
                        {
                        y_i = (int(m) - int(j)) * (int(m) - int(j));
                        for (n = 0; n < side; n++)
                            {
                            x_i = (int(n) - int(k)) * (int(n) - int(k));
                            index = pre_index + (l * side) + m;
                            distance = sqrt(x_i + y_i + z_i);
                            if (distance != 0)
                                reverse_relative_lenghts[index] = distance;
                            }
                        }
                    }
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Wavefunction_square(T* wavefunction_1, T* probabilities, unsigned int lenght_order)
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
T basis_set_calculations<T>::Probabilities_lenght(T* probabilities, unsigned int lenght_order, int x, int y, int z)
    { // Calculate average electron distance from point (x, y, z) for potential calculation
    unsigned int side, size;
    unsigned int i, j, k;
    T lenght;
    T l;
    T m;
    T n;
    T rsqrt;
    
    lenght = 0.00;
    side = (2 * lenght_order + 1);
    size = side * side * side;
    j = 0;
    k = 0;
    
    for (i = 0; i < side; i++)
        {
        l = (i - lenght_order + z) * (i - lenght_order + z);
        for (j = 0; j < side; j++)
            {
            m = (j - lenght_order + y) * (j - lenght_order + y);
            for (k = 0; k < side; k++)
                {
                n = (k - lenght_order + x) * (k - lenght_order + x);
                if (l != 0 or m != 0 or n != 0)
                    {
                    rsqrt = 1/sqrt(l + m + n);
                    lenght = lenght + (rsqrt * T(lenght_order)/T(vector_lenght) * probabilities[(i * side * side) +
                    (j * side) + k]);
                    }
                }
            }
        }
    lenght = 1/lenght * Hartree_lenght; 
    return(lenght);
    }
template <typename T>
T basis_set_calculations<T>::Probabilities_lenght(T* probabilities, unsigned int lenght_order)
    { // Calculate electron distance from null point for potential calculation
    unsigned int side;
    unsigned int size;
    unsigned int i, j, k;
    unsigned int pre_index;
    T lenght;
    T rsqrt;
    T* distances = results.lenghts[0];
    
    lenght = 0.00;
    side = (2 * lenght_order + 1);
    size = side * side * side;
    j = 0;
    k = 0;
    
    for (i = 0; i < side; i++)
        {
        for (j = 0; j < side; j++)
            {
            for (k = 0; k < side; k++)
                {
                pre_index = i * side * side + j * side;
                if (distances[pre_index + k] != 0)
                    {
                    rsqrt = 1/distances[pre_index + k];
                    lenght = lenght + (rsqrt  *
                    probabilities[(i * side * side) + (j * side) + k]);
                    }
                }
            }
        }
    lenght = 1/lenght * T(vector_lenght)/T(lenght_order) * Hartree_lenght;
    return(lenght);
    }
template <typename T>
T basis_set_calculations<T>::Probabilities_thread(T* Probabilities, unsigned int lenght_order, T* lenght)
    { // Calculate electron distance from null point for potential calculation for multithread calls
    lenght[0] = Probabilities_lenght(Probabilities, lenght_order);
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Gradient_thread(T* Gradient_1, T* wavefunction_2, unsigned int lenght_order)
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
T basis_set_calculations<T>::Integral_overlap(T* wavefunction_1, T* wavefunction_2, T* result,
unsigned int lenght_order, T x, T y, T z) 
    { // Calculate the absolute value of overleap integral for multiplying the energy differences from fock matrix
    unsigned int side;
    unsigned int size;
    unsigned int i;
    unsigned int x_contraction, y_contraction, z_contraction;
    int j;
    T* overlap_array;
    T overlap[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    T overlap_sum;
    
    side = (2 * lenght_order + 1);
    x_contraction = abs(x * lenght_order/vector_lenght);
    y_contraction = abs(y * lenght_order/vector_lenght);
    z_contraction = abs(z * lenght_order/vector_lenght);
    if (x_contraction >= side or y_contraction >= side  or z_contraction >= side)
        return(-1);
    
    size = (side - x_contraction) * (side - y_contraction) * (side - z_contraction);
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
        
    if (x == 0 and y == 0 and z == 0)
        Wavefunction_multiply(wavefunction_1, wavefunction_2, overlap_array, lenght_order);
    else
        Wavefunction_multiply(wavefunction_1, wavefunction_2, overlap_array, lenght_order, x, y, z);
    // vectorisation code    
    for (i = 0; i + 7 < size; i = i + 8) // Multiplying of Wavefunction_1 and Wavefunction_1
        {
        overlap[0] = overlap[0] + overlap_array[i];
        overlap[1] = overlap[1] + overlap_array[i + 1];
        overlap[2] = overlap[2] + overlap_array[i + 2];
        overlap[3] = overlap[3] + overlap_array[i + 3];
        overlap[4] = overlap[4] + overlap_array[i + 4];
        overlap[5] = overlap[5] + overlap_array[i + 5];
        overlap[6] = overlap[6] + overlap_array[i + 6];
        overlap[7] = overlap[7] + overlap_array[i + 7]; 
        }
    for (i = size - (size % 8); i < (size); i++)
        overlap[0] = overlap[0] + overlap_array[i];
    
    overlap_sum = overlap[0] + overlap[1] + overlap[2] + overlap[3] + overlap[4] + overlap[5] + overlap[6] + overlap[7];
    // end of vectorisation code    
    delete[] overlap_array;
    if ((not (isnan(overlap_sum))) and (not (isinf(overlap_sum)))) // Check for NaN and inf values
        result[0] = overlap_sum;
    
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Integrate_Integral_overlap(T* wavefunction_1, T* wavefunction_2, T* result,
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
inline T basis_set_calculations<T>::Integral_coulombic(T radius_1, T radius_2, T distance,  T* result, bool spin_bonded)
    { // Calculate the coulombic integral - potential energy of electron to proton of other atoms, aproximate
    T constant;
    T radius;
    T efective_lenght;
    
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
    efective_lenght = sqrt((radius * radius) + (distance * distance));
    if ((not (isnan(constant/efective_lenght))) and (not (isinf(constant/efective_lenght)))) // Check for NaN and inf values
        result[0] = constant/efective_lenght;
    
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Integrate_Integral_coulombic(T* density_1, T* density_2, T* result,
unsigned int lenght_order, T x, T y, T z) 
    { // Calculate the absolute value of overleap integral for multiplying the energy differences from fock matrix
    unsigned int side;
    unsigned int i, j, k, l, m, n;
    unsigned int i_range, j_range, k_range, l_range, m_range, n_range;
    unsigned int i_min, j_min, k_min, l_min, m_min, n_min;
    unsigned int i_max, j_max, k_max, l_max, m_max, n_max;
    unsigned int x_contraction, y_contraction, z_contraction, x_side, y_side, z_side;
    unsigned int pre_index, pre_index_2;
    int x_shift, y_shift, z_shift;
    int x_2_y_2;
    
    T constant;
    T radius;
    T point_value_distance;
    T point_value;
    T average_lenght;
    T* reverse_relative_lenghts = nullptr;
    
    if (small_results.relative_lenghts.size() > 0 and small_results.lenght_orders.size() > 0)
        if (small_results.lenght_orders[0] == lenght_order)
            reverse_relative_lenghts = small_results.relative_lenghts[0];
    
    side = (2 * lenght_order + 1);
    x_shift = x * lenght_order/vector_lenght;
    y_shift = y * lenght_order/vector_lenght;
    z_shift = z * lenght_order/vector_lenght;
    x_contraction = abs(x_shift);
    y_contraction = abs(y_shift);
    z_contraction = abs(z_shift);
    x_side = side;
    y_side = side;
    z_side = side;
    
    // detecting a effective x, y and z ranges of integration
    i_range = 0;
    j_range = 0;
    k_range = 0;
    l_range = 0;
    m_range = 0;
    n_range = 0;
    for (k = 0; k < z_side; k++)
            for (j = 0; j < y_side; j++)
                for (i = 0; i < x_side; i++)
                    {
                    if (density_1[i + j * x_side + k * x_side * y_side] > 0)
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
    for (n = 0; n < z_side; n++)
            for (m = 0; m < y_side; m++)
                for (l = 0; l < x_side; l++)
                    {
                    if (density_2[l + m * x_side + n * x_side * y_side] > 0)
                        {
                        if (l > lenght_order)
                            if (l > l_range + lenght_order)
                                l_range = l - lenght_order;
                        else
                            if (lenght_order > l_range + l)
                                l_range = lenght_order - l;
                        if (m > lenght_order)
                            if (m > m_range + lenght_order)
                                m_range = m - lenght_order;
                        else
                            if (lenght_order > m_range + m)
                                m_range = lenght_order - m;
                        if (n > lenght_order)
                            if (n > n_range + lenght_order)
                                n_range = n - lenght_order;
                        else
                            if (lenght_order > n_range + n)
                                n_range = lenght_order - n;
                        }
                    }
    i_max = lenght_order + i_range;
    j_max = lenght_order + j_range;
    k_max = lenght_order + k_range;
    l_max = lenght_order + l_range;
    m_max = lenght_order + m_range;
    n_max = lenght_order + n_range;
    i_min = lenght_order - i_range;
    j_min = lenght_order - j_range;
    k_min = lenght_order - k_range;
    l_min = lenght_order - l_range;
    m_min = lenght_order - m_range;
    n_min = lenght_order - n_range;
    
    // integration
    average_lenght = 0;
    if (x == 0 and y == 0 and z == 0 and reverse_relative_lenghts != nullptr)
        {
        // optimized integration with zero of coordinate distance
        for (k = k_min; k < k_max; k++)
            for (j = j_min; j < j_max; j++)
                for (i = i_min; i < i_max; i++)
                    for (n = n_min; n < n_max; n++)
                        {
                        pre_index = k * lenght_order * lenght_order * lenght_order * lenght_order * lenght_order +
                        j * lenght_order * lenght_order * lenght_order * lenght_order +
                        i * lenght_order * lenght_order * lenght_order + n * lenght_order * lenght_order;
                        for (m = m_min; m < m_max; m++)
                            {
                            pre_index_2 = m * lenght_order;
                            for (l = l_min; l < l_max; l++)
                                if (i != (l + x_shift) or j != (m + y_shift) or k != (n + z_shift))
                                    {
                                    point_value = reverse_relative_lenghts[pre_index + pre_index_2 + l] *
                                    density_1[i + j * x_side + k * x_side * y_side] *
                                    density_2[l + m * x_side + n * x_side * y_side];
                                    average_lenght = average_lenght + point_value;
                                    }
                            }
                        }
        average_lenght = 1/average_lenght * T(lenght_order)/T(vector_lenght) * Hartree_lenght;
        }
    else
        { // standard integration with non-zero of coordinate distance
        for (k = k_min; k < k_max; k++)
            for (j = j_min; j < j_max; j++)
                for (i = i_min; i < i_max; i++)
                    for (n = n_min; n < n_max; n++)
                        {
                        x_2_y_2 = (int(j) - int(m) - int(y_shift)) * (int(j) - int(m) - int(y_shift)) +
                        (int(k) - int(n) - int(z_shift)) * (int(k) - int(n) - int(z_shift));
                        for (m = m_min; m < m_max; m++)
                            for (l = l_min; l < l_max; l++)
                                if (i != (l + x_shift) or j != (m + y_shift) or k != (n + z_shift))
                                    {
                                    point_value_distance = sqrt((int(i) - int(l) - int(x_shift)) * (int(i) - int(l) - int(x_shift)) +
                                    x_2_y_2);
                                    point_value = 1/(point_value_distance * vector_lenght/lenght_order) *
                                    density_1[i + j * x_side + k * x_side * y_side] *
                                    density_2[l + m * x_side + n * x_side * y_side];
                                    average_lenght = average_lenght + point_value;
                                    }
                            }
        average_lenght = 1/average_lenght  * Hartree_lenght;
        }
    
    constant = e*e/(4*Pi*E0);
    radius = average_lenght;
    
    if ((not (isnan(constant/radius))) and (not (isinf(constant/radius)))) // Check for NaN and inf values
        result[0] = result[0] + constant/radius;
    
    return(average_lenght);
    }
template <typename T>
T basis_set_calculations<T>::Integral_nucleus_atraction(T probabilities_lenght, T multiplier,  T* result, T* lenght, unsigned int Z)
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
T basis_set_calculations<T>::Integrate_Integral_nucleus_atraction(T* probabilities,
T* result, T* lenght, unsigned int lenght_order, T lenght_x, T lenght_y, T lenght_z, unsigned int Z)
    { // Calculate the nucleus attraction integral - potential energy of electron to atom nucleus
    int lenght_X, lenght_Y, lenght_Z; // Z is equal the number of protons - number of eletrons without coulombic integral
    T constant;
    T radius;
    
    constant = -e * e/(4 * Pi * E0);
    lenght_X = lenght_x * lenght_order/vector_lenght;
    lenght_Y = lenght_y * lenght_order/vector_lenght;
    lenght_Z = lenght_z * lenght_order/vector_lenght; 
    radius = Probabilities_lenght(probabilities, lenght_order, lenght_X, lenght_Y, lenght_Z);
    lenght[0] = radius;
        
    if ((not (isnan(constant * Z/radius))) and (not (isinf(constant * Z/radius))))
        result[0] = result[0] + constant * Z/radius; // Check for NaN and inf values
    
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Integral_kinetic(T* Gradient_1, T* Gradient_2, T* result, unsigned int lenght_order,
T d_x, T d_y, T d_z)
    {
    unsigned int i, j, k;
    unsigned int side;
    unsigned int x_side, y_side, z_side;
    unsigned int x_contraction, y_contraction, z_contraction;
    unsigned int x_1_min, y_1_min, z_1_min;
    unsigned int x_2_min, y_2_min, z_2_min;
    unsigned int x_condition, x_condition_2;
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
    if (d_x != 0 or d_y != 0 or d_z != 0)
        for (i = 0; i < (side - z_contraction); i++)
            {
            for (j = 0; j < (side - y_contraction); j++)
                {
                for (k = 0; k + 7 < x_condition; k += 8)
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
                for (k = x_condition_2; k < x_condition; k++)
                    {
                    result_array[0] = result_array[0] + 
                    Gradient_1[((i + z_1_min) * side * side) + ((j + y_1_min) * side) + (k + x_1_min)] *
                    Gradient_2[((i + z_2_min) * side * side) + ((j + y_2_min) * side) + (k + x_2_min)];
                    }
                }
            }
    else
        for (i = 0; i < (side); i++)
            {
            for (j = 0; j < (side); j++)
                {
                for (k = 0; k + 7 < x_condition; k += 8)
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
                for (k = x_condition_2; k < x_condition; k++)
                    {
                    result_array[0] = result_array[0] + Gradient_1[(i * side * side) + (j * side) + k] *
                    Gradient_2[(i * side * side) + (j * side) + k];
                    }
                }
            }
    result[0] = (result_array[0] + result_array[1] + result_array[2] + result_array[3] + result_array[4] + result_array[5]
    + result_array[6] + result_array[7]) * (-h * h)/(8 * Pi * Pi * me);
        // end of vectorisation code
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Rydberg_energy(unsigned int Z, unsigned int n)
    {
    T energy;
    
    energy = -(Z * Z * Hartree_energy_constant/(2 * n * n));
    return(energy);
    }
template <typename T>
T basis_set_calculations<T>::Spin_moment_energy(T s, T B0)
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
T basis_set_calculations<T>::Orbital_magnetic_field(T potential_energy, T radius, int l)
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
T basis_set_calculations<T>::Quantum_numbers_to_orbitals(unsigned int n, unsigned int l, int fulness,
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
T basis_set_calculations<T>::Atoms_to_valence_orbitals(string atom, atom_orbitals* atom_orbitals_PTR)
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
T basis_set_calculations<T>::Create_atomic_wavefunctions(atom_orbitals *atom_orbitals_PTR, atom_wavefunctions *atom_wavefunctions,
    unsigned int size_order,T x,T y,T z)
    { // create or modify atom_wavefunctions list of wavefunctions, probabilities and efective lenghts according to lenght multipliers
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
    atom_wavefunctions->wavefunction_lenght_multipliers.reserve(count_electrons);
    atom_wavefunctions->wavefunction_coefficients.reserve(count_electrons);
    atom_wavefunctions->efective_radius_base.reserve(count_electrons);
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
T basis_set_calculations<T>::Create_bond_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions_1,
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
    unsigned int indexes_1[max_electrons], indexes_2[max_electrons];
    unsigned int indexes_3[max_electrons], indexes_4[max_electrons];
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
T basis_set_calculations<T>::Sum_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions_1, atom_wavefunctions *atom_wavefunctions_2)
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
T basis_set_calculations<T>::Create_nuclear_atraction_integral_matrix(T* matrix, T* nucleuses, unsigned int order,
atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j, k;
    unsigned int count_orbitals;
    unsigned int sum_electrons;
    unsigned int matrix_shift;
    unsigned int ind_nuc_size;
    unsigned int index[max_electrons]; // Index of electrons positions forcomputing nuclear atraction integrals
    unsigned int ind_nuc[max_electrons + 1]; // Index of nucleuses
    
    T* probabilities[max_electrons];
    unsigned int lenght_orders[max_electrons];
    T wavefunction_coefficients[max_electrons];
    T wavefunction_lenght_multipliers[max_electrons];
    T efective_radius_base[max_electrons];
    unsigned int count_electrons[max_electrons];
    unsigned int Z[max_electrons];
    T x[max_electrons];
    T y[max_electrons];
    T z[max_electrons];
    T n[max_electrons];
    int bonding[max_electrons];
    int spin_paired[max_electrons];
    T spins[max_electrons];
    
    thread t61, t62, t63, t64, t65, t66, t67, t68, t69, t70, t71, t72, t73, t74, t75;
    bool restriction;
    bool t69_flag, t70_flag, t71_flag, t72_flag, t73_flag, t74_flag, t75_flag;
    
    t69_flag = false;
    t70_flag = false;
    t71_flag = false;
    t72_flag = false;
    t73_flag = false;
    t74_flag = false;
    t75_flag = false;
    restriction = true;
    matrix_shift = 2;
    
    sum_electrons = atom_wavefunctions->n.size();
    count_orbitals = 0;
    
    if (sum_electrons > max_electrons)
        return(-1);
    
    for (i = 0; i < sum_electrons; i++)
        {
        probabilities[i] = atom_wavefunctions->probabilities[i];
        lenght_orders[i] = atom_wavefunctions->lenght_orders[i];
        wavefunction_coefficients[i] = atom_wavefunctions->wavefunction_coefficients[i];
        wavefunction_lenght_multipliers[i] = atom_wavefunctions->wavefunction_lenght_multipliers[i];
        efective_radius_base[i] = atom_wavefunctions->efective_radius_base[i];
        count_electrons[i] = atom_wavefunctions->count_electrons[i];
        Z[i] = atom_wavefunctions->reduced_Z[i];
        x[i] = atom_wavefunctions->x[i];
        y[i] = atom_wavefunctions->y[i];
        z[i] = atom_wavefunctions->z[i];
        n[i] = atom_wavefunctions->n[i];
        spins[i] = atom_wavefunctions->spins[i];
        spin_paired[i] = atom_wavefunctions->spin_paired[i];
        bonding[i] = atom_wavefunctions->bonding[i];
        }
    for (i = 0; i < (order * order); i++) // initialisation of matrix
        {
        matrix[i] = 0;
        nucleuses[i] = 0;
        }
    // closed-shell basis set method optimalization code
    for (i = 0; i < sum_electrons ; i++) // for restricted and unrestricted basis set method
        {
        if (spin_paired[i] == -1)
            {
            restriction = false;
            matrix_shift = 1;
            }
        }
    for (i = 0; i < sum_electrons ; i++)
        if (wavefunction_coefficients[i] != 0)
            {
            if (i == 0)
                {
                index[count_orbitals] = i;
                count_orbitals++;
                }
            if (i > 0)
                {
                if ((spins[i] == 0.5) or (bonding[i] >= 0) or restriction == false)
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
        for (j = 0; (j + 7) < count_orbitals; j = j + 8)
            {
            if ((x[index[j]] - x[ind_nuc[i]] != 0) or (y[index[j]] - y[ind_nuc[i]] != 0) or
            (z[index[j]] - z[ind_nuc[i]] != 0) and ((n[index[j]]/(Z[index[j]] *
            wavefunction_lenght_multipliers[index[j]])) > (vector_lenght/lenght_orders[index[j]])))
                t61 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this, probabilities[index[j]],
                matrix + (index[j] * (1 + order)), nucleuses_distances + (index[j] * order + ind_nuc[i]),lenght_orders[index[j]],
                x[index[j]] - x[ind_nuc[i]], y[index[j]] - y[ind_nuc[i]], z[index[j]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
            else
                Integral_nucleus_atraction(efective_radius_base[index[j]], wavefunction_lenght_multipliers[index[j]],
                matrix + (index[j] * (1 + order)), nucleuses_distances + (index[j] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            
            if (x[index[j + 1]] - x[ind_nuc[i]] != 0 or y[index[j + 1]] - y[ind_nuc[i]] != 0 or
            z[index[j + 1]] - z[ind_nuc[i]] != 0 and ((n[index[j + 1]]/(Z[index[j + 1]] *
            wavefunction_lenght_multipliers[index[j + 1]])) > (vector_lenght/lenght_orders[index[j + 1]])))
                t62 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this, probabilities[index[j + 1]],
                matrix + (index[j + 1] * (1 + order)), nucleuses_distances + (index[j + 1] * order + ind_nuc[i]),
                lenght_orders[index[j + 1]], x[index[j + 1]] - x[ind_nuc[i]], y[index[j + 1]] - y[ind_nuc[i]],
                z[index[j + 1]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
            else
                Integral_nucleus_atraction(efective_radius_base[index[j + 1]], wavefunction_lenght_multipliers[index[j + 1]],
                matrix + (index[j + 1] * (1 + order)), nucleuses_distances + (index[j + 1] * order + ind_nuc[i]), Z[ind_nuc[i]]);
                
            if (x[index[j + 2]] - x[ind_nuc[i]] != 0 or y[index[j + 2]] - y[ind_nuc[i]] != 0 or
            z[index[j + 2]] - z[ind_nuc[i]] != 0 and ((n[index[j + 2]]/(Z[index[j + 2]] *
            wavefunction_lenght_multipliers[index[j + 2]])) > (vector_lenght/lenght_orders[index[j + 2]])))
                t63 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this, probabilities[index[j + 2]],
                matrix + (index[j + 2] * (1 + order)), nucleuses_distances + (index[j + 2] * order + ind_nuc[i]),
                lenght_orders[index[j + 2]], x[index[j + 2]] - x[ind_nuc[i]], y[index[j + 2]] - y[ind_nuc[i]],
                z[index[j + 2]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
            else
                Integral_nucleus_atraction(efective_radius_base[index[j + 2]], wavefunction_lenght_multipliers[index[j + 2]],
                matrix + (index[j + 2] * (1 + order)), nucleuses_distances + (index[j + 2] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            
            if (x[index[j + 3]] - x[ind_nuc[i]] != 0 or y[index[j + 3]] - y[ind_nuc[i]] != 0 or
            z[index[j + 3]] - z[ind_nuc[i]] != 0 and ((n[index[j + 3]]/(Z[index[j + 3]] *
            wavefunction_lenght_multipliers[index[j + 3]])) > (vector_lenght/lenght_orders[index[j + 3]])))
                t64 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this, probabilities[index[j + 3]],
                matrix + (index[j + 3] * (1 + order)), nucleuses_distances + (index[j + 3] * order + ind_nuc[i]),
                lenght_orders[index[j + 3]], x[index[j + 3]] - x[ind_nuc[i]], y[index[j + 3]] - y[ind_nuc[i]],
                z[index[j + 3]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
            else
                Integral_nucleus_atraction(efective_radius_base[index[j + 3]], wavefunction_lenght_multipliers[index[j + 3]],
                matrix + (index[j + 3] * (1 + order)), nucleuses_distances + (index[j + 3] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            
            if (x[index[j + 4]] - x[ind_nuc[i]] != 0 or y[index[j + 4]] - y[ind_nuc[i]] != 0 or
            z[index[j + 4]] - z[ind_nuc[i]] != 0 and ((n[index[j + 4]]/(Z[index[j + 4]] *
            wavefunction_lenght_multipliers[index[j + 4]])) > (vector_lenght/lenght_orders[index[j + 4]])))
                t65 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this, probabilities[index[j + 4]],
                matrix + (index[j + 4] * (1 + order)), nucleuses_distances + (index[j + 4] * order + ind_nuc[i]),
                lenght_orders[index[j + 4]], x[index[j + 4]] - x[ind_nuc[i]], y[index[j + 4]] - y[ind_nuc[i]],
                z[index[j + 4]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
            else
                Integral_nucleus_atraction(efective_radius_base[index[j + 4]], wavefunction_lenght_multipliers[index[j + 4]],
                matrix + (index[j + 4] * (1 + order)), nucleuses_distances + (index[j + 4] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            
            if (x[index[j + 5]] - x[ind_nuc[i]] != 0 or y[index[j + 5]] - y[ind_nuc[i]] != 0 or
            z[index[j + 5]] - z[ind_nuc[i]] != 0 and ((n[index[j + 5]]/(Z[index[j + 5]] *
            wavefunction_lenght_multipliers[index[j + 5]])) > (vector_lenght/lenght_orders[index[j + 5]])))
                t66 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this, probabilities[index[j + 5]],
                matrix + (index[j + 5] * (1 + order)), nucleuses_distances + (index[j + 5] * order + ind_nuc[i]),
                lenght_orders[index[j + 5]], x[index[j + 5]] - x[ind_nuc[i]], y[index[j + 5]] - y[ind_nuc[i]],
                z[index[j + 5]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
            else
                Integral_nucleus_atraction(efective_radius_base[index[j + 5]], wavefunction_lenght_multipliers[index[j + 5]],
                matrix + (index[j + 5] * (1 + order)), nucleuses_distances + (index[j + 5] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            
            if (x[index[j + 6]] - x[ind_nuc[i]] != 0 or y[index[j + 6]] - y[ind_nuc[i]] != 0 or
            z[index[j + 6]] - z[ind_nuc[i]] != 0 and ((n[index[j + 6]]/(Z[index[j + 6]] *
            wavefunction_lenght_multipliers[index[j + 6]])) > (vector_lenght/lenght_orders[index[j + 6]])))
                t67 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this, probabilities[index[j + 6]],
                matrix + (index[j + 6] * (1 + order)), nucleuses_distances + (index[j + 6] * order + ind_nuc[i]),
                lenght_orders[index[j + 6]], x[index[j + 6]] - x[ind_nuc[i]], y[index[j + 6]] - y[ind_nuc[i]],
                z[index[j + 6]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
            else
                Integral_nucleus_atraction(efective_radius_base[index[j + 6]], wavefunction_lenght_multipliers[index[j + 6]],
                matrix + (index[j + 6] * (1 + order)), nucleuses_distances + (index[j + 6] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            
            if (x[index[j + 7]] - x[ind_nuc[i]] != 0 or y[index[j + 7]] - y[ind_nuc[i]] != 0 or
            z[index[j + 7]] - z[ind_nuc[i]] != 0 and ((n[index[j + 7]]/(Z[index[j + 7]] *
            wavefunction_lenght_multipliers[index[j + 7]])) > (vector_lenght/lenght_orders[index[j + 7]])))
                t68 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this, probabilities[index[j + 7]],
                matrix + (index[j + 7] * (1 + order)), nucleuses_distances + (index[j + 7] * order + ind_nuc[i]),
                lenght_orders[index[j + 7]], x[index[j + 7]] - x[ind_nuc[i]], y[index[j + 7]] - y[ind_nuc[i]],
                z[index[j + 7]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
            else
                Integral_nucleus_atraction(efective_radius_base[index[j + 7]], wavefunction_lenght_multipliers[index[j + 7]],
                matrix + (index[j + 7] * (1 + order)), nucleuses_distances + (index[j + 6] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            
            if (t61.joinable())
                t61.join();
            if (t62.joinable())
                t62.join();
            if (t63.joinable())
                t63.join();
            if (t64.joinable())
                t64.join();
            if (t65.joinable())
                t65.join();
            if (t66.joinable())
                t66.join();
            if (t67.joinable())
                t67.join();
            if (t68.joinable())
                t68.join();
            }
        if (count_orbitals % 8 >= 7)
            {
            if (x[index[count_orbitals - 7]] - x[ind_nuc[i]] != 0 or y[index[count_orbitals - 7]] - y[ind_nuc[i]] != 0
            or z[index[count_orbitals - 7]] - z[ind_nuc[i]] != 0  and ((n[index[count_orbitals - 7]]
            /(Z[index[count_orbitals - 7]] * wavefunction_lenght_multipliers[index[count_orbitals - 7]])) >
            (vector_lenght/lenght_orders[index[count_orbitals - 7]]))) {
                t69 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this,
                probabilities[index[count_orbitals - 7]], matrix + (index[count_orbitals - 7] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 7] * order + ind_nuc[i]), lenght_orders[index[count_orbitals - 7]],
                x[index[count_orbitals - 7]] - x[ind_nuc[i]], y[index[count_orbitals - 7]] - y[ind_nuc[i]],
                z[index[count_orbitals - 7]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
                t69_flag = true;
                }
            else
                Integral_nucleus_atraction(efective_radius_base[index[count_orbitals - 7]],
                wavefunction_lenght_multipliers[index[count_orbitals - 7]],  matrix + (index[count_orbitals - 7] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 7] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            }
        if (count_orbitals % 8 >= 6)
            {
            if (x[index[count_orbitals - 6]] - x[ind_nuc[i]] != 0 or y[index[count_orbitals - 6]] - y[ind_nuc[i]] != 0
            or z[index[count_orbitals - 6]] - z[ind_nuc[i]] != 0  and ((n[index[count_orbitals - 6]]
            /(Z[index[count_orbitals - 6]] * wavefunction_lenght_multipliers[index[count_orbitals - 6]])) >
            (vector_lenght/lenght_orders[index[count_orbitals - 6]]))) {
                t70 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this,
                probabilities[index[count_orbitals - 6]], matrix + (index[count_orbitals - 6] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 6] * order + ind_nuc[i]), lenght_orders[index[count_orbitals - 6]],
                x[index[count_orbitals - 6]] - x[ind_nuc[i]], y[index[count_orbitals - 6]] - y[ind_nuc[i]],
                z[index[count_orbitals - 6]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
                t70_flag = true;
                }
            else
                Integral_nucleus_atraction(efective_radius_base[index[count_orbitals - 6]],
                wavefunction_lenght_multipliers[index[count_orbitals - 6]],  matrix + (index[count_orbitals - 6] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 6] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            }
        if (count_orbitals % 8 >= 5)
            {
            if (x[index[count_orbitals - 5]] - x[ind_nuc[i]] != 0 or y[index[count_orbitals - 5]] - y[ind_nuc[i]] != 0
            or z[index[count_orbitals - 5]] - z[ind_nuc[i]] != 0  and ((n[index[count_orbitals - 5]]
            /(Z[index[count_orbitals - 5]] * wavefunction_lenght_multipliers[index[count_orbitals - 5]])) >
            (vector_lenght/lenght_orders[index[count_orbitals - 5]]))) {
                t71 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this,
                probabilities[index[count_orbitals - 5]], matrix + (index[count_orbitals - 5] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 5] * order + ind_nuc[i]), lenght_orders[index[count_orbitals - 5]],
                x[index[count_orbitals - 5]] - x[ind_nuc[i]], y[index[count_orbitals - 5]] - y[ind_nuc[i]],
                z[index[count_orbitals - 5]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
                t71_flag = true;
                }
            else
                Integral_nucleus_atraction(efective_radius_base[index[count_orbitals - 5]],
                wavefunction_lenght_multipliers[index[count_orbitals - 5]],  matrix + (index[count_orbitals - 5] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 5] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            }
        if (count_orbitals % 8 >= 4)
            {
            if (x[index[count_orbitals - 4]] - x[ind_nuc[i]] != 0 or y[index[count_orbitals - 4]] - y[ind_nuc[i]] != 0
            or z[index[count_orbitals - 4]] - z[ind_nuc[i]] != 0  and ((n[index[count_orbitals - 4]]
            /(Z[index[count_orbitals - 4]] * wavefunction_lenght_multipliers[index[count_orbitals - 4]])) >
            (vector_lenght/lenght_orders[index[count_orbitals - 4]]))) {
                t72 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this,
                probabilities[index[count_orbitals - 4]], matrix + (index[count_orbitals - 4] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 4] * order + ind_nuc[i]), lenght_orders[index[count_orbitals - 4]],
                x[index[count_orbitals - 4]] - x[ind_nuc[i]], y[index[count_orbitals - 4]] - y[ind_nuc[i]],
                z[index[count_orbitals - 4]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
                t72_flag = true;
                }
            else
                Integral_nucleus_atraction(efective_radius_base[index[count_orbitals - 4]],
                wavefunction_lenght_multipliers[index[count_orbitals - 4]],  matrix + (index[count_orbitals - 4] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 4] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            }
        if (count_orbitals % 8 >= 3)
            {
            if (x[index[count_orbitals - 3]] - x[ind_nuc[i]] != 0 or y[index[count_orbitals - 3]] - y[ind_nuc[i]] != 0
            or z[index[count_orbitals - 3]] - z[ind_nuc[i]] != 0  and ((n[index[count_orbitals - 3]]
            /(Z[index[count_orbitals - 3]] * wavefunction_lenght_multipliers[index[count_orbitals - 3]])) >
            (vector_lenght/lenght_orders[index[count_orbitals - 3]]))) {
                t73 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this,
                probabilities[index[count_orbitals - 3]], matrix + (index[count_orbitals - 3] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 3] * order + ind_nuc[i]), lenght_orders[index[count_orbitals - 3]],
                x[index[count_orbitals - 3]] - x[ind_nuc[i]], y[index[count_orbitals - 3]] - y[ind_nuc[i]],
                z[index[count_orbitals - 3]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
                t73_flag = true;
                }
            else
                Integral_nucleus_atraction(efective_radius_base[index[count_orbitals - 3]],
                wavefunction_lenght_multipliers[index[count_orbitals - 3]],  matrix + (index[count_orbitals - 3] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 3] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            }
        if (count_orbitals % 8 >= 2)
            {
            if (x[index[count_orbitals - 2]] - x[ind_nuc[i]] != 0 or y[index[count_orbitals - 2]] - y[ind_nuc[i]] != 0
            or z[index[count_orbitals - 2]] - z[ind_nuc[i]] != 0  and ((n[index[count_orbitals - 2]]
            /(Z[index[count_orbitals - 2]] * wavefunction_lenght_multipliers[index[count_orbitals - 2]])) >
            (vector_lenght/lenght_orders[index[count_orbitals - 2]]))) {
                t74 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this,
                probabilities[index[count_orbitals - 2]], matrix + (index[count_orbitals - 2] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 2] * order + ind_nuc[i]), lenght_orders[index[count_orbitals - 2]],
                x[index[count_orbitals - 2]] - x[ind_nuc[i]], y[index[count_orbitals - 2]] - y[ind_nuc[i]],
                z[index[count_orbitals - 2]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
                t74_flag = true;
                }
            else
                Integral_nucleus_atraction(efective_radius_base[index[count_orbitals - 2]],
                wavefunction_lenght_multipliers[index[count_orbitals - 2]],  matrix + (index[count_orbitals - 2] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 2] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            }
        if (count_orbitals % 8 >= 1)
            {
            if (x[index[count_orbitals - 1]] - x[ind_nuc[i]] != 0 or y[index[count_orbitals - 1]] - y[ind_nuc[i]] != 0
            or z[index[count_orbitals - 1]] - z[ind_nuc[i]] != 0 and ((n[index[count_orbitals - 1]]
            /(Z[index[count_orbitals - 1]] * wavefunction_lenght_multipliers[index[count_orbitals - 1]])) >
            (vector_lenght/lenght_orders[index[count_orbitals - 1]]))) {
                t75 = thread(&basis_set_calculations::Integrate_Integral_nucleus_atraction, this,
                probabilities[index[count_orbitals - 1]], matrix + (index[count_orbitals - 1] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 1] * order + ind_nuc[i]), lenght_orders[index[count_orbitals - 1]],
                x[index[count_orbitals - 1]] - x[ind_nuc[i]], y[index[count_orbitals - 1]] - y[ind_nuc[i]],
                z[index[count_orbitals - 1]] - z[ind_nuc[i]], Z[ind_nuc[i]]);
                t75_flag = true;
                }
            else
                Integral_nucleus_atraction(efective_radius_base[index[count_orbitals - 1]],
                wavefunction_lenght_multipliers[index[count_orbitals - 1]],  matrix + (index[count_orbitals - 1] * (1 + order)),
                nucleuses_distances + (index[count_orbitals - 1] * order + ind_nuc[i]), Z[ind_nuc[i]]);
            }
        if (t69_flag == true)
            {
            t69.join();
            t69_flag = false;
            }
        if (t70_flag == true)
            {
            t70.join();
            t70_flag = false;
            }
        if (t71_flag == true)
            {
            t71.join();
            t71_flag = false;
            }
        if (t72_flag == true)
            {
            t72.join();
            t72_flag = false;
            }
        if (t73_flag == true)
            {
            t73.join();
            t73_flag = false;
            }
        if (t74_flag == true)
            {
            t74.join();
            t74_flag = false;
            }
        if (t75_flag == true)
            {
            t75.join();
            t75_flag = false;
            }
        // end of multithreading code
        for(j = 0; j < order; j++) // copying values into colums of nucleuses
            {
            nucleuses[i + (j * order)] = matrix[j * (order + 1)] * wavefunction_coefficients[j]
            * wavefunction_coefficients[j]; // copy diagonal of matrix
            
            for (k = 0; k < i; k++)
                nucleuses[i + (j * order)] = nucleuses[i + (j * order)] - nucleuses[k + (j * order)];
                
            } // subtract previous values in the row
        }
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
T basis_set_calculations<T>::Create_coulombic_integral_matrix(T* matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions,
small_atom_wavefunctions *small_atom_wavefunctions)
    {
    unsigned int i, j, k;
    unsigned int count_electrons, count_orbitals;
    unsigned int index_size, index_2_size;
    unsigned int small_atom_wavefunctions_size;
    unsigned int index[max_electrons];
    unsigned int index_2_array[max_electrons];
    
    T wavefunction_coefficients[max_electrons];
    T efective_radius_base[max_electrons];
    T wavefunction_lenght_multipliers[max_electrons];
    T lenght_orders[max_electrons];
    T bonding[max_electrons];
    T n[max_electrons];
    T l[max_electrons];
    T spins[max_electrons];
    int spin_paired[max_electrons];
    T x[max_electrons];
    T y[max_electrons];
    T z[max_electrons];
    T* small_probabilities[max_electrons];
    T* small_rel_lenghts;
    T distance, radius_1, radius_2;
    bool restriction;
    bool spin_bonded;
    thread t76, t77, t78, t79, t80, t81, t82, t83, t84, t85, t86, t87, t88, t89, t90;
    bool t84_flag, t85_flag, t86_flag, t87_flag, t88_flag, t89_flag, t90_flag;
    
    t84_flag = false;
    t85_flag = false;
    t86_flag = false;
    t87_flag = false;
    t88_flag = false;
    t89_flag = false;
    t90_flag = false;
    count_electrons = atom_wavefunctions->n.size();
    small_atom_wavefunctions_size = small_atom_wavefunctions->n.size();
    count_orbitals = 0;
    restriction = true;
    
    for (i = 0; i < count_electrons ; i++)
        {
        wavefunction_coefficients[i] = atom_wavefunctions->wavefunction_coefficients[i];
        efective_radius_base[i] = atom_wavefunctions->efective_radius_base[i];
        wavefunction_lenght_multipliers[i] = atom_wavefunctions->wavefunction_lenght_multipliers[i];
        if (i < small_atom_wavefunctions_size)
            lenght_orders[i] = small_atom_wavefunctions->lenght_orders[i];
        
        bonding[i] = atom_wavefunctions->bonding[i];
        n[i] = atom_wavefunctions->n[i];
        l[i] = atom_wavefunctions->l[i];
        spins[i] = atom_wavefunctions->spins[i];
        spin_paired[i] = atom_wavefunctions->spin_paired[i];
        x[i] = atom_wavefunctions->x[i];
        y[i] = atom_wavefunctions->y[i];
        z[i] = atom_wavefunctions->z[i];
        if (i < small_atom_wavefunctions_size)
            small_probabilities[i] = small_atom_wavefunctions->probabilities[i];
        }
    if (small_atom_wavefunctions_size > 0)
        small_rel_lenghts = small_atom_wavefunctions->relative_lenghts[0];
    for (i = 0; i < (order * order); i++) // Initializing matrix array
        matrix[i] = 0;
        
    for (i = 0; i < count_electrons; i++) // Using regression curve for s1 - s1 integrals
        {
        for (j = i + 1; j < count_electrons; j++)
            if ((n[i] == 1 and n[j] == 1) or ((n[i] == 1 or n[j] == 1)
            and (spin_paired[i] >= 0 or spin_paired[j] >= 0) and (x[j] - x[i] == 0 and y[j] - y[i] == 0 and z[j] - z[i] == 0))
            or (spin_paired[i] == j and (x[j] - x[i] == 0 and y[j] - y[i] == 0 and z[j] - z[i] == 0)))
                {
                radius_1 = efective_radius_base[i] * wavefunction_lenght_multipliers[i];
                radius_2 = efective_radius_base[j] * wavefunction_lenght_multipliers[j];
                distance = sqrt(((x[j] - x[i]) * (x[j] - x[i])) + ((y[j] - y[i]) * (y[j] - y[i])) + ((z[j] - z[i]) * (z[j] - z[i])))
                * Hartree_lenght;
                if ((spin_paired[i] == j or (l[i] == 0 and l[j] == 0)) and (x[j] - x[i] == 0 and y[j] - y[i] == 0 and z[j] - z[i] == 0))
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
        for (j = 0; (j + 7) < index_2_size; j = j + 8)
            { // calculate coulombic integrals multithread
            t76 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[j]], &matrix[index[i] + (index_2_array[j] * order)], lenght_orders[index_2_array[j]],
            x[index_2_array[j]] - x[index[i]], y[index_2_array[j]] - y[index[i]], z[index_2_array[j]] - z[index[i]]);
            t77 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[j + 1]], &matrix[index[i] + (index_2_array[j + 1] * order)],
            lenght_orders[index_2_array[j + 1]], x[index_2_array[j + 1]] - x[index[i]], y[index_2_array[j + 1]] - y[index[i]],
            z[index_2_array[j + 1]] - z[index[i]]);
            t78 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[j + 2]], &matrix[index[i] + (index_2_array[j + 2] * order)],
            lenght_orders[index_2_array[j + 2]], x[index_2_array[j + 2]] - x[index[i]], y[index_2_array[j + 2]] - y[index[i]],
            z[index_2_array[j + 2]] - z[index[i]]);
            t79 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[j + 3]], &matrix[index[i] + (index_2_array[j + 3] * order)],
            lenght_orders[index_2_array[j + 3]], x[index_2_array[j + 3]] - x[index[i]], y[index_2_array[j + 3]] - y[index[i]],
            z[index_2_array[j + 3]] - z[index[i]]);
            t80 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[j + 4]], &matrix[index[i] + (index_2_array[j + 4] * order)],
            lenght_orders[index_2_array[j + 4]], x[index_2_array[j + 4]] - x[index[i]], y[index_2_array[j + 4]] - y[index[i]],
            z[index_2_array[j + 4]] - z[index[i]]);
            t81 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[j + 5]], &matrix[index[i] + (index_2_array[j + 5] * order)],
            lenght_orders[index_2_array[j + 5]], x[index_2_array[j + 5]] - x[index[i]], y[index_2_array[j + 5]] - y[index[i]],
            z[index_2_array[j + 5]] - z[index[i]]);
            t82 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[j + 6]], &matrix[index[i] + (index_2_array[j + 6] * order)],
            lenght_orders[index_2_array[j + 6]], x[index_2_array[j + 6]] - x[index[i]], y[index_2_array[j + 6]] - y[index[i]],
            z[index_2_array[j + 6]] - z[index[i]]);
            t83 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[j + 7]], &matrix[index[i] + (index_2_array[j + 7] * order)],
            lenght_orders[index_2_array[j + 7]], x[index_2_array[j + 7]] - x[index[i]], y[index_2_array[j + 7]] - y[index[i]],
            z[index_2_array[j + 7]] - z[index[i]]);
            
            t76.join();
            t77.join();
            t78.join();
            t79.join();
            t80.join();
            t81.join();
            t82.join();
            t83.join();
            }
    
        if (index_2_size % 8 >= 7)
            {
            t84 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[index_2_size - 7]], &matrix[index[i] + (index_2_array[index_2_size - 7] * order)],
            lenght_orders[index_2_array[index_2_size - 7]], x[index_2_array[index_2_size - 7]] - x[index[i]],
            y[index_2_array[index_2_size - 7]] - y[index[i]], z[index_2_array[index_2_size - 7]] - z[index[i]]);
            t84_flag = true;
            }
        if (index_2_size % 8 >= 6)
            {
            t85 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[index_2_size - 6]], &matrix[index[i] + (index_2_array[index_2_size - 6] * order)],
            lenght_orders[index_2_array[index_2_size - 6]], x[index_2_array[index_2_size - 6]] - x[index[i]],
            y[index_2_array[index_2_size - 6]] - y[index[i]], z[index_2_array[index_2_size - 6]] - z[index[i]]);
            t85_flag = true;
            }
        if (index_2_size % 8 >= 5)
            {
            t86 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[index_2_size - 5]], &matrix[index[i] + (index_2_array[index_2_size - 5] * order)],
            lenght_orders[index_2_array[index_2_size - 5]], x[index_2_array[index_2_size - 5]] - x[index[i]],
            y[index_2_array[index_2_size - 5]] - y[index[i]], z[index_2_array[index_2_size - 5]] - z[index[i]]);
            t86_flag = true;
            }
        if (index_2_size % 8 >= 4)
            {
            t87 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[index_2_size - 4]], &matrix[index[i] + (index_2_array[index_2_size - 4] * order)],
            lenght_orders[index_2_array[index_2_size - 4]], x[index_2_array[index_2_size - 4]] - x[index[i]],
            y[index_2_array[index_2_size - 4]] - y[index[i]], z[index_2_array[index_2_size - 4]] - z[index[i]]);
            t87_flag = true;
            }
        if (index_2_size % 8 >= 3)
            {
            t88 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[index_2_size - 3]], &matrix[index[i] + (index_2_array[index_2_size - 3] * order)],
            lenght_orders[index_2_array[index_2_size - 3]], x[index_2_array[index_2_size - 3]] - x[index[i]],
            y[index_2_array[index_2_size - 3]] - y[index[i]], z[index_2_array[index_2_size - 3]] - z[index[i]]);
            t88_flag = true;
            }
        if (index_2_size % 8 >= 2)
            {
            t89 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[index_2_size - 2]], &matrix[index[i] + (index_2_array[index_2_size - 2] * order)],
            lenght_orders[index_2_array[index_2_size - 2]], x[index_2_array[index_2_size - 2]] - x[index[i]],
            y[index_2_array[index_2_size - 2]] - y[index[i]], z[index_2_array[index_2_size - 2]] - z[index[i]]);
            t89_flag = true;
            }
        if (index_2_size % 8 >= 1)
            {
            t90 = thread(&basis_set_calculations::Integrate_Integral_coulombic, this, small_probabilities[index[i]],
            small_probabilities[index_2_array[index_2_size - 1]], &matrix[index[i] + (index_2_array[index_2_size - 1] * order)],
            lenght_orders[index_2_array[index_2_size - 1]], x[index_2_array[index_2_size - 1]] - x[index[i]],
            y[index_2_array[index_2_size - 1]] - y[index[i]], z[index_2_array[index_2_size - 1]] - z[index[i]]);
            t90_flag = true;
            }
        if (t84_flag == true)
            {
            t84.join();
            t84_flag = false;
            }
        if (t85_flag == true)
            {
            t85.join();
            t85_flag = false;
            }
        if (t86_flag == true)
            {
            t86.join();
            t86_flag = false;
            }
        if (t87_flag == true)
            {
            t87.join();
            t87_flag = false;
            }
        if (t88_flag == true)
            {
            t88.join();
            t88_flag = false;
            }
        if (t89_flag == true)
            {
            t89.join();
            t89_flag = false;
            }
        if (t90_flag == true)
            {
            t90.join();
            t90_flag = false;
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
                matrix[j + (i * order)] = matrix[j + (i * order)] * wavefunction_coefficients[i] * wavefunction_coefficients[i]
                * wavefunction_coefficients[j] * wavefunction_coefficients[j];
            else
                { // for orbitals with wavefunction_coefficients > 1 are coulombic and exchange energies cancelled
                if (wavefunction_coefficients[j] > 1)
                    matrix[j + (i * order)] = matrix[j + (i * order)] * wavefunction_coefficients[i] * wavefunction_coefficients[i];
                if (wavefunction_coefficients[i] > 1)
                    matrix[j + (i * order)] = matrix[j + (i * order)] * wavefunction_coefficients[j] * wavefunction_coefficients[j];
                }
            }
    for (i = 0; i < order; i++) // copy calculated efective_lenghts upper diagonal
        for (j = 0; j < i; j++)
            matrix[(j * order) + i] = matrix[(i * order) + j];
            
    for (i = 0; i < order; i++) // copying 0 to diagonal
        matrix[i * (order + 1)] = 0;
        
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Create_overlap_integral_matrix(T* matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j, k;
    unsigned int count_electrons;
    unsigned int count_orbitals;
    unsigned int index_size;
    unsigned int index_2_size;
    
    unsigned int index[max_electrons];
    unsigned int index_2_array[max_electrons];
    T* wavefunctions[max_electrons];
    unsigned int lenght_orders[max_electrons];
    T wavefunction_coefficients[max_electrons];
    T x[max_electrons];
    T y[max_electrons];
    T z[max_electrons];
    T spins[max_electrons];
    unsigned int l[max_electrons];
    int m[max_electrons];
    int spin_paired[max_electrons];
    int bonding[max_electrons];
    bool t99_flag, t100_flag, t101_flag, t102_flag, t103_flag, t104_flag, t105_flag;
    bool restriction;
    thread t91, t92, t93, t94, t95, t96, t97, t98, t99, t100, t101, t102, t103, t104, t105;
    
    t99_flag = false;
    t100_flag = false;
    t101_flag = false;
    t102_flag = false;
    t103_flag = false;
    t104_flag = false;
    t105_flag = false;
    restriction = true;
    count_electrons = atom_wavefunctions->n.size();
    count_orbitals = 0;
    
    for (i = 0; i < count_electrons ; i++)
        {
        wavefunctions[i] = atom_wavefunctions->wavefunctions[i];
        lenght_orders[i] = atom_wavefunctions->lenght_orders[i];
        wavefunction_coefficients[i] = atom_wavefunctions->wavefunction_coefficients[i];
        x[i] = atom_wavefunctions->x[i];
        y[i] = atom_wavefunctions->y[i];
        z[i] = atom_wavefunctions->z[i];
        spins[i] = atom_wavefunctions->spins[i];
        l[i] = atom_wavefunctions->l[i];
        m[i] = atom_wavefunctions->m[i];
        spin_paired[i] = atom_wavefunctions->spin_paired[i];
        bonding[i] = atom_wavefunctions->bonding[i];
        }
    for (i = 0; i < (order * order); i++) // initializing matrix array
        matrix[i] = 0;
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
        for (j = 0; (j + 7) < index_2_size; j = j + 8)
            { // calculate overlap integrals multithread
            t91 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[j]], &matrix[index[i] + (index_2_array[j] * order)], lenght_orders[index_2_array[j]],
            x[index_2_array[j]] - x[index[i]], y[index_2_array[j]] - y[index[i]], z[index_2_array[j]] - z[index[i]]);
            t92 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[j + 1]], &matrix[index[i] + (index_2_array[j + 1] * order)],
            lenght_orders[index_2_array[j + 1]], x[index_2_array[j + 1]] - x[index[i]], y[index_2_array[j + 1]] - y[index[i]],
            z[index_2_array[j + 1]] - z[index[i]]);
            t93 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[j + 2]], &matrix[index[i] + (index_2_array[j + 2] * order)],
            lenght_orders[index_2_array[j + 2]], x[index_2_array[j + 2]] - x[index[i]], y[index_2_array[j + 2]] - y[index[i]],
            z[index_2_array[j + 2]] - z[index[i]]);
            t94 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[j + 3]], &matrix[index[i] + (index_2_array[j + 3] * order)],
            lenght_orders[index_2_array[j + 3]], x[index_2_array[j + 3]] - x[index[i]], y[index_2_array[j + 3]] - y[index[i]],
            z[index_2_array[j + 3]] - z[index[i]]);
            t95 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[j + 4]], &matrix[index[i] + (index_2_array[j + 4] * order)],
            lenght_orders[index_2_array[j + 4]], x[index_2_array[j + 4]] - x[index[i]], y[index_2_array[j + 4]] - y[index[i]],
            z[index_2_array[j + 4]] - z[index[i]]);
            t96 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[j + 5]], &matrix[index[i] + (index_2_array[j + 5] * order)],
            lenght_orders[index_2_array[j + 5]], x[index_2_array[j + 5]] - x[index[i]], y[index_2_array[j + 5]] - y[index[i]],
            z[index_2_array[j + 5]] - z[index[i]]);
            t97 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[j + 6]], &matrix[index[i] + (index_2_array[j + 6] * order)],
            lenght_orders[index_2_array[j + 6]], x[index_2_array[j + 6]] - x[index[i]], y[index_2_array[j + 6]] - y[index[i]],
            z[index_2_array[j + 6]] - z[index[i]]);
            t98 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[j + 7]], &matrix[index[i] + (index_2_array[j + 7] * order)],
            lenght_orders[index_2_array[j + 7]], x[index_2_array[j + 7]] - x[index[i]], y[index_2_array[j + 7]] - y[index[i]],
            z[index_2_array[j + 7]] - z[index[i]]);
            
            t91.join();
            t92.join();
            t93.join();
            t94.join();
            t95.join();
            t96.join();
            t97.join();
            t98.join();
            }
        if (index_2_size % 8 >= 7)
            {
            t99 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[index_2_size - 7]], &matrix[index[i] + (index_2_array[index_2_size - 7] * order)],
            lenght_orders[index_2_array[index_2_size - 7]], x[index_2_array[index_2_size - 7]] - x[index[i]],
            y[index_2_array[index_2_size - 7]] - y[index[i]], z[index_2_array[index_2_size - 7]] - z[index[i]]);
            t99_flag = true;
            }
        if (index_2_size % 8 >= 6)
            {
            t100 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[index_2_size - 6]], &matrix[index[i] + (index_2_array[index_2_size - 6] * order)],
            lenght_orders[index_2_array[index_2_size - 6]], x[index_2_array[index_2_size - 6]] - x[index[i]],
            y[index_2_array[index_2_size - 6]] - y[index[i]], z[index_2_array[index_2_size - 6]] - z[index[i]]);
            t100_flag = true;
            }
        if (index_2_size % 8 >= 5)
            {
            t101 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[index_2_size - 5]], &matrix[index[i] + (index_2_array[index_2_size - 5] * order)],
            lenght_orders[index_2_array[index_2_size - 5]], x[index_2_array[index_2_size - 5]] - x[index[i]],
            y[index_2_array[index_2_size - 5]] - y[index[i]], z[index_2_array[index_2_size - 5]] - z[index[i]]);
            t101_flag = true;
            }
        if (index_2_size % 8 >= 4)
            {
            t102 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[index_2_size - 4]], &matrix[index[i] + (index_2_array[index_2_size - 4] * order)],
            lenght_orders[index_2_array[index_2_size - 4]], x[index_2_array[index_2_size - 4]] - x[index[i]],
            y[index_2_array[index_2_size - 4]] - y[index[i]], z[index_2_array[index_2_size - 4]] - z[index[i]]);
            t102_flag = true;
            }
        if (index_2_size % 8 >= 3)
            {
            t103 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[index_2_size - 3]], &matrix[index[i] + (index_2_array[index_2_size - 3] * order)],
            lenght_orders[index_2_array[index_2_size - 3]], x[index_2_array[index_2_size - 3]] - x[index[i]],
            y[index_2_array[index_2_size - 3]] - y[index[i]], z[index_2_array[index_2_size - 3]] - z[index[i]]);
            t103_flag = true;
            }
        if (index_2_size % 8 >= 2)
            {
            t104 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[index_2_size - 2]], &matrix[index[i] + (index_2_array[index_2_size - 2] * order)],
            lenght_orders[index_2_array[index_2_size - 2]], x[index_2_array[index_2_size - 2]] - x[index[i]],
            y[index_2_array[index_2_size - 2]] - y[index[i]], z[index_2_array[index_2_size - 2]] - z[index[i]]);
            t104_flag = true;
            }
        if (index_2_size % 8 >= 1)
            {
            t105 = thread(&basis_set_calculations::Integral_overlap, this, wavefunctions[index[i]],
            wavefunctions[index_2_array[index_2_size - 1]], &matrix[index[i] + (index_2_array[index_2_size - 1] * order)],
            lenght_orders[index_2_array[index_2_size - 1]], x[index_2_array[index_2_size - 1]] - x[index[i]],
            y[index_2_array[index_2_size - 1]] - y[index[i]], z[index_2_array[index_2_size - 1]] - z[index[i]]);
            t105_flag = true;
            }
        if (t99_flag == true)
            {
            t99.join();
            t99_flag = false;
            }
        if (t100_flag == true)
            {
            t100.join();
            t100_flag = false;
            }
        if (t101_flag == true)
            {
            t101.join();
            t101_flag = false;
            }
        if (t102_flag == true)
            {
            t102.join();
            t102_flag = false;
            }
        if (t103_flag == true)
            {
            t103.join();
            t103_flag = false;
            }
        if (t104_flag == true)
            {
            t104.join();
            t104_flag = false;
            }
        if (t105_flag == true)
            {
            t105.join();
            t105_flag = false;
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
            
    for (i = 0; i < order; i++) // copy calculated efective_lenghts upper diagonal
        for (j = 0; j < i; j++)
            matrix[(j * order) + i] = matrix[(i * order) + j];
            
    for (i = 0; i < order; i++) // copying 0 to diagonal
        matrix[i * (order + 1)] = 0;

    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Calculate_resonance_integral_matrix(T* overlap_matrix, T* overlap_efective_lenght_integral_matrix,
T* resonance_integral_matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions,
small_atom_wavefunctions *small_atom_wavefunctions)
    {
    unsigned int i, j, k;
    unsigned int count_electrons, count_small_wavefunctions;
    unsigned int count_interactions;
    
    T spins[max_electrons];
    T x[max_electrons];
    T y[max_electrons];
    T z[max_electrons];
    T n[max_electrons];
    unsigned int electron_numbers[max_electrons];
    T efective_radius_base[max_electrons];
    T wavefunction_lenght_multipliers[max_electrons];
    T wavefunction_coefficients[max_electrons];
    int bonding[max_electrons];
    int antibonding[max_electrons];
    unsigned int  small_electron_numbers[max_electrons];
    unsigned int  small_lenght_orders[max_electrons];
    T* small_wavefunctions[max_electrons];
    T small_x[max_electrons];
    T small_y[max_electrons];
    T small_z[max_electrons];
    unsigned int index[max_electrons];
    unsigned int index_size;
    T radius, radius_1, radius_2;
    T efective_lenght;
    T constant;
    T overlap;
    T pow_overrlaps;
    
    thread t106, t107, t108, t109, t110, t111, t112, t113, t114, t115, t116, t117, t118, t119, t120;
    bool t114_flag, t115_flag, t116_flag, t117_flag, t118_flag, t119_flag, t120_flag;
    
    
    t114_flag = false;
    t115_flag = false;
    t116_flag = false;
    t117_flag = false;
    t118_flag = false;
    t119_flag = false;
    t120_flag = false;
    
    count_electrons = atom_wavefunctions->n.size();
    count_small_wavefunctions = small_atom_wavefunctions->electron_numbers.size();
    for (i = 0; i < count_electrons; i++)
        {
        spins[i] = atom_wavefunctions->spins[i];
        x[i] = atom_wavefunctions->x[i];
        y[i] = atom_wavefunctions->y[i];
        z[i] = atom_wavefunctions->z[i];
        n[i] = atom_wavefunctions->n[i];
        electron_numbers[i] = atom_wavefunctions->electron_numbers[i];
        efective_radius_base[i] = atom_wavefunctions->efective_radius_base[i];
        wavefunction_lenght_multipliers[i] = atom_wavefunctions->wavefunction_lenght_multipliers[i];
        wavefunction_coefficients[i] = atom_wavefunctions->wavefunction_coefficients[i];
        bonding[i] = atom_wavefunctions->bonding[i];
        antibonding[i] = atom_wavefunctions->antibonding[i];
        }
    for (i = 0; i < count_small_wavefunctions; i++)
        {
        small_electron_numbers[i] = small_atom_wavefunctions->electron_numbers[i];
        small_lenght_orders[i] = small_atom_wavefunctions->lenght_orders[i];
        small_wavefunctions[i] = small_atom_wavefunctions->wavefunctions[i];
        small_x[i] = small_atom_wavefunctions->x[i];
        small_y[i] = small_atom_wavefunctions->y[i];
        small_z[i] = small_atom_wavefunctions->z[i];
        }
    constant = e*e/(4*Pi*E0);
    
    for (i = 0; i < order * order; i++)
        overlap_efective_lenght_integral_matrix[i] = 0;
    for (i = 0; i < order * order; i++)
        resonance_integral_matrix[i] = 0;
    
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
                    radius_1 = efective_radius_base[i]/wavefunction_lenght_multipliers[i];
                    radius_2 = efective_radius_base[j]/wavefunction_lenght_multipliers[j];
                    if (radius_1 >= radius_2)
                        radius = radius_1 + (Phi + 1)/3.00 * radius_2;
                    else
                        radius = radius_2 + (Phi + 1)/3.00 * radius_1;  
                    efective_lenght = radius * (1.481178577 - (0.700616171 * abs(overlap_matrix[(i * order) + j]))
                    + (0.224068020 * abs(overlap_matrix[(i * order) + j]) * abs(overlap_matrix[(i * order) + j])));
                    overlap_efective_lenght_integral_matrix[(i * order) +j] = efective_lenght;
                    }
                else
                    { // Integrating of efective lenght of wavefunction overlap for resonance integral from small wavefunction set
                    index[k] = j;
                    k++;
                    }
                }
            }
        index_size = k;
        // multithreading code  
        for (j = 0; (j + 7) < index_size; j = j + 8)
            { // calculate overlap integrals multithread
            t106 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[j]], &overlap_efective_lenght_integral_matrix[(i * order) + index[j]], small_lenght_orders[i],
            small_x[index[j]] - small_x[i], small_y[index[j]] - small_y[i], small_z[index[j]] - small_z[i]);
            t107 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[j + 1]], &overlap_efective_lenght_integral_matrix[(i * order) + index[j + 1]],
            small_lenght_orders[i], small_x[index[j + 1]] - small_x[i], small_y[index[j] + 1] - small_y[i],
            small_z[index[j] + 1] - small_z[i]);
            t108 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[j + 2]], &overlap_efective_lenght_integral_matrix[(i * order) + index[j + 2]],
            small_lenght_orders[i], small_x[index[j + 2]] - small_x[i], small_y[index[j] + 2] - small_y[i],
            small_z[index[j] + 2] - small_z[i]);
            t109 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[j + 3]], &overlap_efective_lenght_integral_matrix[(i * order) + index[j + 3]],
            small_lenght_orders[i], small_x[index[j + 3]] - small_x[i], small_y[index[j] + 3] - small_y[i],
            small_z[index[j] + 3] - small_z[i]);
            t110 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[j + 4]], &overlap_efective_lenght_integral_matrix[(i * order) + index[j + 4]],
            small_lenght_orders[i], small_x[index[j + 4]] - small_x[i], small_y[index[j] + 4] - small_y[i],
            small_z[index[j] + 4] - small_z[i]);
            t111 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[j + 5]], &overlap_efective_lenght_integral_matrix[(i * order) + index[j + 5]],
            small_lenght_orders[i], small_x[index[j + 5]] - small_x[i], small_y[index[j] + 5] - small_y[i],
            small_z[index[j] + 5] - small_z[i]);
            t112 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[j + 6]], &overlap_efective_lenght_integral_matrix[(i * order) + index[j + 6]],
            small_lenght_orders[i], small_x[index[j + 6]] - small_x[i], small_y[index[j] + 6] - small_y[i],
            small_z[index[j] + 6] - small_z[i]);
            t113 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[j + 7]], &overlap_efective_lenght_integral_matrix[(i * order) + index[j + 7]],
            small_lenght_orders[i], small_x[index[j + 7]] - small_x[i], small_y[index[j] + 7] - small_y[i],
            small_z[index[j] + 7] - small_z[i]);
            
            t106.join();
            t107.join();
            t108.join();
            t109.join();
            t110.join();
            t111.join();
            t112.join();
            t113.join();
            }
        if (index_size % 8 >= 7)
            {
            t114 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[index_size - 7]], &overlap_efective_lenght_integral_matrix[(i * order) + index[index_size - 7]],
            small_lenght_orders[i], small_x[index[index_size - 7]] - small_x[i], small_y[index[index_size - 7]] - small_y[i],
            small_z[index[index_size - 7]] - small_z[i]);
            t114_flag = true;
            }
        if (index_size % 8 >= 6)
            {
            t115 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[index_size - 6]], &overlap_efective_lenght_integral_matrix[(i * order) + index[index_size - 6]],
            small_lenght_orders[i], small_x[index[index_size - 6]] - small_x[i], small_y[index[index_size - 6]] - small_y[i],
            small_z[index[index_size - 6]] - small_z[i]);
            t115_flag = true;
            }
        if (index_size % 8 >= 5)
            {
            t116 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[index_size - 5]], &overlap_efective_lenght_integral_matrix[(i * order) + index[index_size - 5]],
            small_lenght_orders[i], small_x[index[index_size - 5]] - small_x[i], small_y[index[index_size - 5]] - small_y[i],
            small_z[index[index_size - 5]] - small_z[i]);
            t116_flag = true;
            }
        if (index_size % 8 >= 4)
            {
            t117 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[index_size - 4]], &overlap_efective_lenght_integral_matrix[(i * order) + index[index_size - 4]],
            small_lenght_orders[i], small_x[index[index_size - 4]] - small_x[i], small_y[index[index_size - 4]] - small_y[i],
            small_z[index[index_size - 4]] - small_z[i]);
            t117_flag = true;
            }
        if (index_size % 8 >= 3)
            {
            t118 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[index_size - 3]], &overlap_efective_lenght_integral_matrix[(i * order) + index[index_size - 3]],
            small_lenght_orders[i], small_x[index[index_size - 3]] - small_x[i], small_y[index[index_size - 3]] - small_y[i],
            small_z[index[index_size - 3]] - small_z[i]);
            t118_flag = true;
            }
        if (index_size % 8 >= 2)
            {
            t119 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[index_size - 2]], &overlap_efective_lenght_integral_matrix[(i * order) + index[index_size - 2]],
            small_lenght_orders[i], small_x[index[index_size - 2]] - small_x[i], small_y[index[index_size - 2]] - small_y[i],
            small_z[index[index_size - 2]] - small_z[i]);
            t119_flag = true;
            }
        if (index_size % 8 >= 1)
            {
            t120 = thread(&basis_set_calculations::Integrate_Integral_overlap, this, small_wavefunctions[i],
            small_wavefunctions[index[index_size - 1]], &overlap_efective_lenght_integral_matrix[(i * order) + index[index_size - 1]],
            small_lenght_orders[i], small_x[index[index_size - 1]] - small_x[i], small_y[index[index_size - 1]] - small_y[i],
            small_z[index[index_size - 1]] - small_z[i]);
            t120_flag = true;
            }
        if (t114_flag == true)
            {
            t114.join();
            t114_flag = false;
            }
        if (t115_flag == true)
            {
            t115.join();
            t115_flag = false;
            }
        if (t116_flag == true)
            {
            t116.join();
            t116_flag = false;
            }
        if (t117_flag == true)
            {
            t117.join();
            t117_flag = false;
            }
        if (t118_flag == true)
            {
            t118.join();
            t118_flag = false;
            }
        if (t119_flag == true)
            {
            t119.join();
            t119_flag = false;
            }
        if (t120_flag == true)
            {
            t120.join();
            t120_flag = false;
            }
        } // end of multithreading code
    for (i = 0; i < order; i++) // copy calculated efective_lenghts upper diagonal
        for (j = 0; j < i; j++)
            overlap_efective_lenght_integral_matrix[(j * order) + i] = overlap_efective_lenght_integral_matrix[(i * order) + j];
        
    for (i = 0; i < order; i++)
        for (j = 0; j < i; j++)
            if ((overlap_matrix[(i * order) + j] != 0) and (bonding[i] >= 0 or antibonding[i] >= 0) and
            (bonding[j] >= 0 or antibonding[j] >= 0) and (x[i] != x[j] or y[i] != y[j] or z[i] != z[j]))
                {
                resonance_integral_matrix[(i * order) + j] = constant/overlap_efective_lenght_integral_matrix[(i * order) + j] *
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
T basis_set_calculations<T>::Calculate_kinetic_integral_matrix(T* matrix, unsigned int order, atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j, k;
    unsigned int count_electrons, count_orbitals;
    unsigned int index_size, index_2_size;
    unsigned int index[max_electrons];
    unsigned int index_2_array[max_electrons];
    int spin_paired[max_electrons];
    
    T* Gradients[max_electrons];
    T wavefunction_coefficients[max_electrons];
    T efective_radius_base[max_electrons];
    T wavefunction_lenght_multipliers[max_electrons];
    T lenght_orders[max_electrons];
    T bonding[max_electrons];
    T n[max_electrons];
    T l[max_electrons];
    T spins[max_electrons];
    T x[max_electrons];
    T y[max_electrons];
    T z[max_electrons];
    T Z[max_electrons];
    T pixel_lenght;
    
    pixel_lenght = T(vector_lenght)/T(order) * Hartree_lenght;
    
    T distance, radius_1, radius_2;
    bool restriction;
    
    thread t121, t122, t123, t124, t125, t126, t127, t128, t129, t130, t131, t132, t133, t134, t135;
    bool t129_flag, t130_flag, t131_flag, t132_flag, t133_flag, t134_flag, t135_flag;
    
    t129_flag = false;
    t130_flag = false;
    t131_flag = false;
    t132_flag = false;
    t133_flag = false;
    t134_flag = false;
    t135_flag = false;
    count_electrons = atom_wavefunctions->n.size();
    count_orbitals = 0;
    restriction = true;
    
    for (i = 0; i < count_electrons ; i++)
        {
        Gradients[i] = atom_wavefunctions->Gradients[i];
        wavefunction_coefficients[i] = atom_wavefunctions->wavefunction_coefficients[i];
        efective_radius_base[i] = atom_wavefunctions->efective_radius_base[i];
        wavefunction_lenght_multipliers[i] = atom_wavefunctions->wavefunction_lenght_multipliers[i];
        lenght_orders[i] = atom_wavefunctions->lenght_orders[i];
        bonding[i] = atom_wavefunctions->bonding[i];
        n[i] = atom_wavefunctions->n[i];
        l[i] = atom_wavefunctions->l[i];
        spins[i] = atom_wavefunctions->spins[i];
        spin_paired[i] = atom_wavefunctions->spin_paired[i];
        x[i] = atom_wavefunctions->x[i];
        y[i] = atom_wavefunctions->y[i];
        z[i] = atom_wavefunctions->z[i];
        Z[i] = atom_wavefunctions->Z[i];
        }
        
    for (i = 0; i < (order * order); i++) // Initializing matrix array
        matrix[i] = 0;

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
        for (j = 0; (j + 7) < index_2_size; j = j + 8)
            { // calculate kinetic integrals multithread
            t121 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[j]], &matrix[index[i] + (index_2_array[j] * order)], lenght_orders[index_2_array[j]],
            x[index_2_array[j]] - x[index[i]], y[index_2_array[j]] - y[index[i]], z[index_2_array[j]] - z[index[i]]);
            t122 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[j + 1]], &matrix[index[i] + (index_2_array[j + 1] * order)],
            lenght_orders[index_2_array[j + 1]], x[index_2_array[j + 1]] - x[index[i]], y[index_2_array[j + 1]] - y[index[i]],
            z[index_2_array[j + 1]] - z[index[i]]);
            t123 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[j + 2]], &matrix[index[i] + (index_2_array[j + 2] * order)],
            lenght_orders[index_2_array[j + 2]], x[index_2_array[j + 2]] - x[index[i]], y[index_2_array[j + 2]] - y[index[i]],
            z[index_2_array[j + 2]] - z[index[i]]);
            t124 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[j + 3]], &matrix[index[i] + (index_2_array[j + 3] * order)],
            lenght_orders[index_2_array[j + 3]], x[index_2_array[j + 3]] - x[index[i]], y[index_2_array[j + 3]] - y[index[i]],
            z[index_2_array[j + 3]] - z[index[i]]);
            t125 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[j + 4]], &matrix[index[i] + (index_2_array[j + 4] * order)],
            lenght_orders[index_2_array[j + 4]], x[index_2_array[j + 4]] - x[index[i]], y[index_2_array[j + 4]] - y[index[i]],
            z[index_2_array[j + 4]] - z[index[i]]);
            t126 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[j + 5]], &matrix[index[i] + (index_2_array[j + 5] * order)],
            lenght_orders[index_2_array[j + 5]], x[index_2_array[j + 5]] - x[index[i]], y[index_2_array[j + 5]] - y[index[i]],
            z[index_2_array[j + 5]] - z[index[i]]);
            t127 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[j + 6]], &matrix[index[i] + (index_2_array[j + 6] * order)],
            lenght_orders[index_2_array[j + 6]], x[index_2_array[j + 6]] - x[index[i]], y[index_2_array[j + 6]] - y[index[i]],
            z[index_2_array[j + 6]] - z[index[i]]);
            t128 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[j + 7]], &matrix[index[i] + (index_2_array[j + 7] * order)],
            lenght_orders[index_2_array[j + 7]], x[index_2_array[j + 7]] - x[index[i]], y[index_2_array[j + 7]] - y[index[i]],
            z[index_2_array[j + 7]] - z[index[i]]);
            
            t121.join();
            t122.join();
            t123.join();
            t124.join();
            t125.join();
            t126.join();
            t127.join();
            t128.join();
            }
        if (index_2_size % 8 >= 7) {
            t129 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[index_2_size - 7]], &matrix[index[i] + (index_2_array[index_2_size - 7] * order)],
            lenght_orders[index_2_array[index_2_size - 7]], x[index_2_array[index_2_size - 7]] - x[index[i]],
            y[index_2_array[index_2_size - 7]] - y[index[i]], z[index_2_array[index_2_size - 7]] - z[index[i]]);
            t129_flag = true;
            }
        if (index_2_size % 8 >= 6) {
            t130 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[index_2_size - 6]], &matrix[index[i] + (index_2_array[index_2_size - 6] * order)],
            lenght_orders[index_2_array[index_2_size - 6]], x[index_2_array[index_2_size - 6]] - x[index[i]],
            y[index_2_array[index_2_size - 6]] - y[index[i]], z[index_2_array[index_2_size - 6]] - z[index[i]]);
            t130_flag = true;
            }
        if (index_2_size % 8 >= 5) {
            t131 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[index_2_size - 5]], &matrix[index[i] + (index_2_array[index_2_size - 5] * order)],
            lenght_orders[index_2_array[index_2_size - 5]], x[index_2_array[index_2_size - 5]] - x[index[i]],
            y[index_2_array[index_2_size - 5]] - y[index[i]], z[index_2_array[index_2_size - 5]] - z[index[i]]);
            t131_flag = true;
            }
        if (index_2_size % 8 >= 4) {
            t132 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[index_2_size - 4]], &matrix[index[i] + (index_2_array[index_2_size - 4] * order)],
            lenght_orders[index_2_array[index_2_size - 4]], x[index_2_array[index_2_size - 4]] - x[index[i]],
            y[index_2_array[index_2_size - 4]] - y[index[i]], z[index_2_array[index_2_size - 4]] - z[index[i]]);
            t132_flag = true;
            }
        if (index_2_size % 8 >= 3) {
            t133 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[index_2_size - 3]], &matrix[index[i] + (index_2_array[index_2_size - 3] * order)],
            lenght_orders[index_2_array[index_2_size - 3]], x[index_2_array[index_2_size - 3]] - x[index[i]],
            y[index_2_array[index_2_size - 3]] - y[index[i]], z[index_2_array[index_2_size - 3]] - z[index[i]]);
            t133_flag = true;
            }
        if (index_2_size % 8 >= 2) {
            t134 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[index_2_size - 2]], &matrix[index[i] + (index_2_array[index_2_size - 2] * order)],
            lenght_orders[index_2_array[index_2_size - 2]], x[index_2_array[index_2_size - 2]] - x[index[i]],
            y[index_2_array[index_2_size - 2]] - y[index[i]], z[index_2_array[index_2_size - 2]] - z[index[i]]);
            t134_flag = true;
            }
        if (index_2_size % 8 >= 1) {
            t135 = thread(&basis_set_calculations::Integral_kinetic, this, Gradients[index[i]],
            Gradients[index_2_array[index_2_size - 1]], &matrix[index[i] + (index_2_array[index_2_size - 1] * order)],
            lenght_orders[index_2_array[index_2_size - 1]], x[index_2_array[index_2_size - 1]] - x[index[i]],
            y[index_2_array[index_2_size - 1]] - y[index[i]], z[index_2_array[index_2_size - 1]] - z[index[i]]);
            t135_flag = true;
            }
        if (t129_flag == true) {
            t129.join();
            t129_flag = false;
            }
        if (t130_flag == true) {
            t130.join();
            t130_flag = false;
            }
        if (t131_flag == true) {
            t131.join();
            t131_flag = false;
            }
        if (t132_flag == true) {
            t132.join();
            t132_flag = false;
            }
        if (t133_flag == true) {
            t133.join();
            t133_flag = false;
            }
        if (t134_flag == true) {
            t134.join();
            t134_flag = false;
            }
        if (t135_flag == true) {
            t135.join();
            t135_flag = false;
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
    for (i = 0; i < order; i++) // copy calculated efective_lenghts upper diagonal
        for (j = 0; j < i; j++)
            matrix[(j * order) + i] = matrix[(i * order) + j];
        
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Calculate_basis_set_matrix(T* nuclear_atraction_integral_matrix, T* coulombic_integral_matrix,
T* resonance_integral_matrix, T* kinetic_integral_matrix, T* basis_set_matrix, 
unsigned int order, atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j;
    
    T spin_orbit_energy;
    T B;
    T distance_nucleuses;
    T radius, radius_1, radius_2;
    T efective_lenght;
    T efective_radius_base[max_electrons];
    T wavefunction_lenght_multipliers[max_electrons];
    T potential_energy[max_electrons];
    
    unsigned int Z[max_electrons];
    unsigned int x[max_electrons];
    unsigned int y[max_electrons];
    unsigned int z[max_electrons];
    unsigned int n[max_electrons];
    unsigned int l[max_electrons];
    int m[max_electrons];
    T spins[max_electrons];
    int spin_paired[max_electrons];
    
    for (i = 0; i < order; i++)
        {
        efective_radius_base[i] = atom_wavefunctions->efective_radius_base[i];
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
                    
                    radius_1 = efective_radius_base[i]/wavefunction_lenght_multipliers[i];
                    radius_2 = efective_radius_base[j]/wavefunction_lenght_multipliers[j];
                    if (radius_1 >= radius_2)
                        radius = radius_1 + (Phi - 1) * radius_2;
                    else
                        radius = radius_2 + (Phi - 1) * radius_1;
    
                    efective_lenght = sqrt((radius * radius) + (distance_nucleuses * distance_nucleuses));
                       
                    B = Orbital_magnetic_field(potential_energy[i], efective_lenght, l[j]);
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
T basis_set_calculations<T>::Calculate_corr_basis_set_matrix(T* basis_set_matrix, T* correction_matrix,
T* corr_basis_set_matrix, unsigned int order)
    {
    unsigned int i;
    
    for (i = 0; i < order * order; i++)
        corr_basis_set_matrix[i] = basis_set_matrix[i] + correction_matrix[i];
    
    return(0);
    }
template <typename T>    
T basis_set_calculations<T>::Calculate_spin_density_matrix(T* overlap_integral_matrix, T* spin_density_matrix, unsigned int order,
atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j;
    T spins[max_electrons];
    
    for (i = 0; i < order; i++)
        spins[i] = atom_wavefunctions->spins[i];
    
    for (i = 0; i < order; i++)
        for (j = 0; j < order; j++)
            spin_density_matrix[(i * order) + j] = abs(overlap_integral_matrix[(i * order) + j]) * spins[i] * 2;
            
    for (i = 0; i < order; i++) // copying values to diagonal
            spin_density_matrix[i * (order + 1)] = spins[i] * 2;
    
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Solve_basis_set_matrix(T* basis_set_matrix, T* overlap_integral_matrix, unsigned int order,
vector<T>* values, atom_wavefunctions *atom_wavefunctions)
    {
    unsigned int i, j, k;
    unsigned int indexes_2[max_electrons];
    
    T alpha, beta, E;
    T multiplier_constant, multiplier;
    T Hamiltonian_parts[max_electrons];
    T beta_sum_rows[max_electrons];
    T overlaps[max_electrons];
    T Eigenvectors[max_electrons];
    T Hamiltonian = 0;
    T new_old_iteration_ratio[] =  {0.50, 0.50,
    0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250,
    0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250,
    0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
    0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
    0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
    0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
    0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,
    0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625,};
    
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
            
        Eigenvectors[i] = E; // determinants into Eigenvectors sorted  according basis_set matrix rows
        }
    for (i = 0; i < order; i++) // calculating Hamiltonian
        if ((not (isnan(Hamiltonian_parts[i])) and (not isinf(Hamiltonian_parts[i])))) // Check for NaN and inf values
            Hamiltonian = Hamiltonian + Hamiltonian_parts[i];
        
    for (i = 0; i < order; i++) // copying to vectors
        {// algoritm of generating new wavefunctions cause doubling of values
        if (atom_wavefunctions->bonding[i] >= 0)
            {
            values->push_back(-abs((Eigenvectors[i] + Eigenvectors[bonding[i]])/2));
            electron_spectra.push_back(abs((Eigenvectors[i] + Eigenvectors[bonding[i]])/2));
            }
        else
            {
            values->push_back(Eigenvectors[i]);
            electron_spectra.push_back(-Eigenvectors[i]);
            }
        }
    for (i = 0; i < order; i++) // variational generation of new wavefunction lenght coefficients for next iteration from roots of matrix
        {
        if ((Z[i] - charge[i]))
            multiplier_constant = 1.00/sqrt(Z[i] - charge[i]);
            
        if (Eigenvectors[i] != 0 and -correction_matrix[i * (1 + order)] != Rydberg_energy(Z[i], n[i]) and constraints[i] == 0)
            {
            multiplier = pow(wavefunction_lenght_multipliers[i], 1 - new_old_iteration_ratio[Z[i] - charge[i]]) *
            pow(sqrt(abs(Rydberg_energy(Z[i], n[i]) /Eigenvectors[i] * multiplier_constant)), new_old_iteration_ratio[Z[i] - charge[i]])
            * Get_relative_Hartree_length(Z[i], n[i]);
            if (multiplier > 1.00/1024 and multiplier < 1024)
                wavefunction_lenght_multipliers[i] = multiplier;
            }
        }
    sort(electron_spectra.begin(),electron_spectra.end());
    return(Hamiltonian);
    }
template <typename T>
T basis_set_calculations<T>::Solve_spin_density_matrix(T* spin_density_matrix, unsigned int order,
    vector<T>* energetic_levels)
    {
    unsigned int i, j;
    unsigned int det_size, spectra_size;
    
    basis_set_Determinant_solver(order, spin_density_matrix);
    det_size = determinants.size();
    spectra_size = spectra_EPR.size();
    energetic_levels->clear();
    
    for (i = 0; i < det_size; i++)
        energetic_levels->push_back(determinants[i]);
    
    return(0);
    }
// end of section 5 - generating matrices of integrals and Fock matrices, section 6: completing basis set method
// and user interface handling
template <typename T>
T basis_set_calculations<T>::Atom_orbitals_generate(string UI_input, atom_orbitals *atom_orbitals_PTR)
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
    Atoms_to_valence_orbitals(input, atom_orbitals_PTR);
    atom_orbitals_PTR->charge = -shift;
    atom_orbitals_PTR->Z = element_number + 1;
    atom_orbitals_PTR->reduced_Z = atom_orbitals_PTR->reduced_Z - shift;
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::Nucleus_repulsive_energy(atom_wavefunctions *atom_wavefunctions)
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
T basis_set_calculations<T>::Generate_atomic_wavefunctions(atom_wavefunctions *atom_wavefunctions,
small_atom_wavefunctions *small_atom_wavefunctions, unsigned int size_order, bool alocate, bool compute_densities)
    { // create or modify atom_wavefunctions list of wavefunctions, probabilities and efective lenghts according to lenght multipliers
    unsigned int i, j, k;
    unsigned int count_orbitals;
    unsigned int count_electrons;
    unsigned int wavefunction_size;
    unsigned int spin_paired_electron_index;
    unsigned int n_i[max_electrons];
    unsigned int l_i[max_electrons];
    int m_i[max_electrons];
    unsigned int Z[max_electrons];
    unsigned int Z_i[max_electrons];
    
    T* pointers_to_lenghts[1];
    T* pointers_to_wavefunctions[max_electrons];
    T* pointers_to_probabilities[max_electrons];
    T* pointers_to_Gradients[max_electrons];
    T* pointer_to_lenght;
    T* pointer_to_wavefunction;
    T* pointer_to_probability;
    T* pointer_to_Gradient;
    T* small_lenghts;
    T* small_relative_lenghts;
    T* small_wavefunctions[max_electrons];
    T* small_probabilities[max_electrons];
    unsigned int small_lenght_order;
    unsigned int small_wavefunction_size;
    unsigned int small_atom_wavefunctions_size;
    unsigned int  small_electron_numbers[max_electrons];
    unsigned int  small_lenght_orders[max_electrons];
    
    T multiplier;
    unsigned int index[max_electrons + 1];
    T multipliers[max_electrons];
    T multiplier_array[max_electrons];
    T efective_radius_base_array[max_electrons];
    T spins[max_electrons];
    int spin_paired[max_electrons];
    int bonding[max_electrons];
    unsigned int constraints[max_electrons];
    unsigned int electron_numbers[max_electrons];
    
    bool restriction;
    bool non_s1_system;
    bool t9_flag, t10_flag, t11_flag, t12_flag, t13_flag, t14_flag, t15_flag;
    bool t24_flag, t25_flag, t26_flag, t27_flag, t28_flag, t29_flag, t30_flag;
    bool t39_flag, t40_flag, t41_flag, t42_flag, t43_flag, t44_flag, t45_flag;
    bool t54_flag, t55_flag, t56_flag, t57_flag, t58_flag, t59_flag, t60_flag;
    thread t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
    thread t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30;
    thread t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42, t43, t44, t45;
    thread t46, t47, t48, t49, t50, t51, t52, t53, t54, t55, t56, t57, t58, t59, t60;
    
    count_orbitals = 0;
    count_electrons = atom_wavefunctions->n.size();
    
    if (count_electrons > max_electrons)
        return(-1);
        
    wavefunction_size = ((2 * size_order) + 1) * ((2 * size_order) + 1) * ((2 * size_order) + 1);
    small_lenght_order = (sqrt(size_order) * 1.5);
    small_wavefunction_size = (2 * small_lenght_order + 1) * (2 * small_lenght_order + 1) * (2 * small_lenght_order + 1);
    t9_flag = false;
    t10_flag = false;
    t11_flag = false;
    t12_flag = false;
    t13_flag = false;
    t14_flag = false;
    t15_flag = false;
    t24_flag = false;
    t25_flag = false;
    t26_flag = false;
    t27_flag = false;
    t28_flag = false;
    t29_flag = false;
    t30_flag = false;
    t39_flag = false;
    t40_flag = false;
    t41_flag = false;
    t42_flag = false;
    t43_flag = false;
    t44_flag = false;
    t45_flag = false;
    t54_flag = false;
    t55_flag = false;
    t56_flag = false;
    t57_flag = false;
    t58_flag = false;
    t59_flag = false;
    t60_flag = false;
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
    if (alocate == true) // generating lists of lenghts for wavefunctions, wavefunctions, probabilities and lenghts multipliers
        {
        for (i = 0; i < count_electrons; i++)
            {
            if (restriction == false)
                {
                if (spins[i] == 0.5)
                    atom_wavefunctions->wavefunction_lenght_multipliers[i] *=  1.005;
                else
                    atom_wavefunctions->wavefunction_lenght_multipliers[i] *=  0.995;
                }
            pointers_to_wavefunctions[i] = nullptr;
            pointers_to_probabilities[i] = nullptr;
            pointers_to_Gradients[i] = nullptr;
            small_wavefunctions[i] = nullptr;
            small_probabilities[i] = nullptr;
            }
        
        pointers_to_lenghts[0] = nullptr;
        small_lenghts = nullptr;
        small_relative_lenghts = nullptr;
        try
            {
            pointer_to_lenght = new T[wavefunction_size];
            pointers_to_lenghts[0] = pointer_to_lenght;
            atom_wavefunctions->lenghts.push_back(pointers_to_lenghts[0]);
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
            for (i = 0; (i < count_orbitals) and (i < 8); i++)
                {
                if (pointers_to_lenghts[i] != nullptr)
                    delete[] pointers_to_lenghts[i];
                }
            if (small_lenghts != nullptr)
                    delete[] small_lenghts;
            
            if (small_relative_lenghts != nullptr)
                    delete[] small_relative_lenghts;
            
            return(-1);
            }
        // Generating lenghts 3D cubes for wavefunction calculations
        Wavefunction_lenghts_generate(pointers_to_lenghts[0], size_order);
        // Generating small_lenghts 3D cubes for wavefunction calculations
        if (non_s1_system == true) {
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
            { // If allocation fail then dealocate alocated memory and return -1;
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
            if (non_s1_system == true)
                for (i = 0; i < count_electrons; i++)
                    {
                    small_electron_numbers[i] = electron_numbers[i];
                    small_lenght_orders[i] = small_lenght_order;
                    small_wavefunctions[i] = new T[small_wavefunction_size];
                    small_probabilities[i] = new T[small_wavefunction_size];
                    small_atom_wavefunctions->electron_numbers.push_back(small_electron_numbers[i]);
                    small_atom_wavefunctions->lenght_orders.push_back(small_lenght_orders[i]);
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
        pointers_to_lenghts[0] = atom_wavefunctions->lenghts[0];
        if (non_s1_system == true)
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
            small_wavefunctions[i] = small_atom_wavefunctions->wavefunctions[i];
            small_probabilities[i] = small_atom_wavefunctions->probabilities[i];
            }
        }
    for (i = 0; i < count_orbitals; i++)
            {
            multipliers[i] = atom_wavefunctions->wavefunction_lenght_multipliers[index[i]];
            multiplier_array[i] = (vector_lenght/size_order) * multipliers[i];
            } // end of closed-shell basis set method optimalization code
    if (non_s1_system == true)
        for (i = 0; i < small_atom_wavefunctions_size; i++) // Generating small wavefunctions
            { 
            Orbitals_to_wavefunctions(small_atom_wavefunctions->n[i], small_atom_wavefunctions->l[i], small_atom_wavefunctions->m[i],
            small_lenght_orders[i], small_wavefunctions[i], small_lenghts,
            small_atom_wavefunctions->Z[i], atom_wavefunctions->wavefunction_lenght_multipliers[i] *
            (vector_lenght/small_lenght_order));
            Wavefunction_square(small_wavefunctions[i], small_probabilities[i], small_lenght_order);
            }
    // multithreading code
    for (i = 0; (i + 7) < count_orbitals; i = i + 8)
        { // multithreading generation of wavefunctions
        t1 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[i], l_i[i], m_i[i], size_order,
        pointers_to_wavefunctions[index[i]], pointers_to_lenghts[0], Z_i[i], multiplier_array[i]);
        t2 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[i + 1], l_i[i + 1], m_i[i + 1], size_order,
        pointers_to_wavefunctions[index[i + 1]], pointers_to_lenghts[0], Z_i[i + 1], multiplier_array[i + 1]);
        t3 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[i + 2], l_i[i + 2], m_i[i + 2], size_order,
        pointers_to_wavefunctions[index[i + 2]], pointers_to_lenghts[0], Z_i[i + 2], multiplier_array[i + 2]);
        t4 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[i + 3], l_i[i + 3], m_i[i + 3], size_order,
        pointers_to_wavefunctions[index[i + 3]], pointers_to_lenghts[0], Z_i[i + 3], multiplier_array[i + 3]);
        t5 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[i + 4], l_i[i + 4], m_i[i + 4], size_order,
        pointers_to_wavefunctions[index[i + 4]], pointers_to_lenghts[0], Z_i[i + 4], multiplier_array[i + 4]);
        t6 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[i + 5], l_i[i + 5], m_i[i + 5], size_order,
        pointers_to_wavefunctions[index[i + 5]], pointers_to_lenghts[0], Z_i[i + 5], multiplier_array[i + 5]);
        t7 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[i + 6], l_i[i + 6], m_i[i + 6], size_order,
        pointers_to_wavefunctions[index[i + 6]], pointers_to_lenghts[0], Z_i[i + 6], multiplier_array[i + 6]);
        t8 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[i + 7], l_i[i + 7], m_i[i + 7], size_order,
        pointers_to_wavefunctions[index[i + 7]], pointers_to_lenghts[0], Z_i[i + 7], multiplier_array[i + 7]);
        t1.join();
        t2.join();
        t3.join();
        t4.join();
        t5.join();
        t6.join();
        t7.join();
        t8.join();
        }
    if (count_orbitals % 8 >= 7) {
        t9 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[count_orbitals - 7],
        l_i[count_orbitals - 7], m_i[count_orbitals - 7], size_order, pointers_to_wavefunctions[index[count_orbitals - 7]],
        pointers_to_lenghts[0], Z_i[count_orbitals - 7], multiplier_array[count_orbitals - 7]);
        t9_flag = true;
        }
    if (count_orbitals % 8 >= 6) {
        t10 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[count_orbitals - 6],
        l_i[count_orbitals - 6], m_i[count_orbitals - 6], size_order, pointers_to_wavefunctions[index[count_orbitals - 6]],
        pointers_to_lenghts[0], Z_i[count_orbitals - 6], multiplier_array[count_orbitals - 6]);
        t10_flag = true;
        }
    if (count_orbitals % 8 >= 5) {
        t11 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[count_orbitals - 5],
        l_i[count_orbitals - 5], m_i[count_orbitals - 5], size_order, pointers_to_wavefunctions[index[count_orbitals - 5]],
        pointers_to_lenghts[0], Z_i[count_orbitals - 5], multiplier_array[count_orbitals - 5]);
        t11_flag = true;
        }
    if (count_orbitals % 8 >= 4) {
        t12 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[count_orbitals - 4],
        l_i[count_orbitals - 4], m_i[count_orbitals - 4], size_order, pointers_to_wavefunctions[index[count_orbitals - 4]],
        pointers_to_lenghts[0], Z_i[count_orbitals - 4], multiplier_array[count_orbitals - 4]);
        t12_flag = true;
        }
    if (count_orbitals % 8 >= 3) {
        t13 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[count_orbitals - 3],
        l_i[count_orbitals - 3], m_i[count_orbitals - 3], size_order, pointers_to_wavefunctions[index[count_orbitals - 3]],
        pointers_to_lenghts[0], Z_i[count_orbitals - 3], multiplier_array[count_orbitals - 3]);
        t13_flag = true;
        }
    if (count_orbitals % 8 >= 2) {
        t14 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[count_orbitals - 2],
        l_i[count_orbitals - 2], m_i[count_orbitals - 2], size_order, pointers_to_wavefunctions[index[count_orbitals - 2]],
        pointers_to_lenghts[0], Z_i[count_orbitals - 2], multiplier_array[count_orbitals - 2]);
        t14_flag = true;
        }
    if (count_orbitals % 8 >= 1) {
        t15 = thread(&basis_set_calculations::Orbitals_to_wavefunctions, this, n_i[count_orbitals - 1],
        l_i[count_orbitals - 1], m_i[count_orbitals - 1], size_order, pointers_to_wavefunctions[index[count_orbitals - 1]],
        pointers_to_lenghts[0], Z_i[count_orbitals - 1], multiplier_array[count_orbitals - 1]);
        t15_flag = true;
        }
    if (t9_flag == true) {
        t9.join();
        t9_flag = false;
        }
    if (t10_flag == true) {
        t10.join();
        t10_flag = false;
        }
    if (t11_flag == true) {
        t11.join();
        t11_flag = false;
        }
    if (t12_flag == true) {
        t12.join();
        t12_flag = false;
        }
    if (t13_flag == true) {
        t13.join();
        t13_flag = false;
        }
    if (t14_flag == true) {
        t14.join();
        t14_flag = false;
        }
    if (t15_flag == true) {
        t15.join();
        t15_flag = false;
        }
    if (compute_densities == true)
        {
        for (i = 0; (i + 7) < count_orbitals; i = i + 8)
            { // multithreading generation of probabilities
            t16 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[i]],
            pointers_to_probabilities[index[i]], size_order);
            t17 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[i + 1]],
            pointers_to_probabilities[index[i + 1]], size_order);
            t18 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[i + 2]],
            pointers_to_probabilities[index[i + 2]], size_order);
            t19 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[i + 3]],
            pointers_to_probabilities[index[i + 3]], size_order);
            t20 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[i + 4]],
            pointers_to_probabilities[index[i + 4]], size_order);
            t21 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[i + 5]],
            pointers_to_probabilities[index[i + 5]], size_order);
            t22 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[i + 6]],
            pointers_to_probabilities[index[i + 6]], size_order);
            t23 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[i + 7]],
            pointers_to_probabilities[index[i + 7]], size_order);
            t16.join();
            t17.join();
            t18.join();
            t19.join();
            t20.join();
            t21.join();
            t22.join();
            t23.join();
            }
        if (count_orbitals % 8 >= 7) {
            t24 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[count_orbitals - 7]],
            pointers_to_probabilities[index[count_orbitals - 7]], size_order);
            t24_flag = true;
            }
        if (count_orbitals % 8 >= 6) {
            t25 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[count_orbitals - 6]],
            pointers_to_probabilities[index[count_orbitals - 6]], size_order);
            t25_flag = true;
            }
        if (count_orbitals % 8 >= 5) {
            t26 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[count_orbitals - 5]],
            pointers_to_probabilities[index[count_orbitals - 5]], size_order);
            t26_flag = true;
            }
        if (count_orbitals % 8 >= 4) {
            t27 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[count_orbitals - 4]],
            pointers_to_probabilities[index[count_orbitals - 4]], size_order);
            t27_flag = true;
            }
        if (count_orbitals % 8 >= 3) {
            t28 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[count_orbitals - 3]],
            pointers_to_probabilities[index[count_orbitals - 3]], size_order);
            t28_flag = true;
            }
        if (count_orbitals % 8 >= 2) {
            t29 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[count_orbitals - 2]],
            pointers_to_probabilities[index[count_orbitals - 2]], size_order);
            t29_flag = true;
            }
        if (count_orbitals % 8 >= 1) {
            t30 = thread(&basis_set_calculations::Wavefunction_square, this, pointers_to_wavefunctions[index[count_orbitals - 1]],
            pointers_to_probabilities[index[count_orbitals - 1]], size_order);
            t30_flag = true;
            }
        if (t24_flag == true) {
            t24.join();
            t24_flag = false;
            }
        if (t25_flag == true) {
            t25.join();
            t25_flag = false;
            }
        if (t26_flag == true) {
            t26.join();
            t26_flag = false;
            }
        if (t27_flag == true) {
            t27.join();
            t27_flag = false;
            }
        if (t28_flag == true) {
            t28.join();
            t28_flag = false;
            }
        if (t29_flag == true) {
            t29.join();
            t29_flag = false;
            }
        if (t30_flag == true) {
            t30.join();
            t30_flag = false;
            }
        }
    if (alocate == true)
        {
        for (i = 0; (i + 7) < count_orbitals; i = i + 8)
            { // multithreading generation of efective_radius_base parameters
            t31 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[i]],
            size_order, efective_radius_base_array + index[i]);
            t32 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[i + 1]],
            size_order, efective_radius_base_array + index[i + 1]);
            t33 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[i + 2]],
            size_order, efective_radius_base_array + index[i + 2]);
            t34 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[i + 3]],
            size_order, efective_radius_base_array + index[i + 3]);
            t35 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[i + 4]],
            size_order, efective_radius_base_array + index[i + 4]);
            t36 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[i + 5]],
            size_order, efective_radius_base_array + index[i + 5]);
            t37 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[i + 6]],
            size_order, efective_radius_base_array + index[i + 6]);
            t38 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[i + 7]],
            size_order, efective_radius_base_array + index[i + 7]);
            t31.join();
            t32.join();
            t33.join();
            t34.join();
            t35.join();
            t36.join();
            t37.join();
            t38.join();
            }
        if (count_orbitals % 8 >= 7) {
            t39 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[count_orbitals - 7]],
            size_order, efective_radius_base_array + index[count_orbitals - 7]);
            t39_flag = true;
            }
        if (count_orbitals % 8 >= 6) {
            t40 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[count_orbitals - 6]],
            size_order, efective_radius_base_array + index[count_orbitals - 6]);
            t40_flag = true;
            }
        if (count_orbitals % 8 >= 5) {
            t41 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[count_orbitals - 5]],
            size_order, efective_radius_base_array + index[count_orbitals - 5]);
            t41_flag = true;
            }
        if (count_orbitals % 8 >= 4) {
            t42 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[count_orbitals - 4]],
            size_order, efective_radius_base_array + index[count_orbitals - 4]);
            t42_flag = true;
            }
        if (count_orbitals % 8 >= 3) {
            t43 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[count_orbitals - 3]],
            size_order, efective_radius_base_array + index[count_orbitals - 3]);
            t43_flag = true;
            }
        if (count_orbitals % 8 >= 2) {
            t44 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[count_orbitals - 2]],
            size_order, efective_radius_base_array + index[count_orbitals - 2]);
            t44_flag = true;
            }
        if (count_orbitals % 8 >= 1) {
            t45 = thread(&basis_set_calculations::Probabilities_thread, this, pointers_to_probabilities[index[count_orbitals - 1]],
            size_order, efective_radius_base_array + index[count_orbitals - 1]);
            t45_flag = true;
            }
        if (t39_flag == true) {
            t39.join();
            t39_flag = false;
            }
        if (t40_flag == true) {
            t40.join();
            t40_flag = false;
            }
        if (t41_flag == true) {
            t41.join();
            t41_flag = false;
            }
        if (t42_flag == true) {
            t42.join();
            t42_flag = false;
            }
        if (t43_flag == true) {
            t43.join();
            t43_flag = false;
            }
        if (t44_flag == true) {
            t44.join();
            t44_flag = false;
            }
        if (t45_flag == true) {
            t45.join();
            t45_flag = false;
            }
        }
    for (i = 0; (i + 7) < count_orbitals; i = i + 8)
        {
        t46 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[i]],
        pointers_to_wavefunctions[index[i]], size_order);
        t47 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[i + 1]],
        pointers_to_wavefunctions[index[i + 1]], size_order);
        t48 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[i + 2]],
        pointers_to_wavefunctions[index[i + 2]], size_order);
        t49 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[i + 3]],
        pointers_to_wavefunctions[index[i + 3]], size_order);
        t50 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[i + 4]],
        pointers_to_wavefunctions[index[i + 4]], size_order);
        t51 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[i + 5]],
        pointers_to_wavefunctions[index[i + 5]], size_order);
        t52 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[i + 6]],
        pointers_to_wavefunctions[index[i + 6]], size_order);
        t53 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[i + 7]],
        pointers_to_wavefunctions[index[i + 7]], size_order);
            
        t46.join();
        t47.join();
        t48.join();
        t49.join();
        t50.join();
        t51.join();
        t52.join();
        t53.join();
        }
    if (count_orbitals % 8 >= 7) {
        t54 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[count_orbitals - 7]],
        pointers_to_wavefunctions[index[count_orbitals - 7]], size_order);
        t54_flag = true;
        }
    if (count_orbitals % 8 >= 6) {
        t55 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[count_orbitals - 6]],
        pointers_to_wavefunctions[index[count_orbitals - 6]], size_order);
        t55_flag = true;
        }
    if (count_orbitals % 8 >= 5) {
        t56 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[count_orbitals - 5]],
        pointers_to_wavefunctions[index[count_orbitals - 5]], size_order);
        t56_flag = true;
        }
    if (count_orbitals % 8 >= 4) {
        t57 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[count_orbitals - 4]],
        pointers_to_wavefunctions[index[count_orbitals - 4]], size_order);
        t57_flag = true;
        }
    if (count_orbitals % 8 >= 3) {
        t58 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[count_orbitals - 3]],
        pointers_to_wavefunctions[index[count_orbitals - 3]], size_order);
        t58_flag = true;
        }
    if (count_orbitals % 8 >= 2) {
        t59 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[count_orbitals - 2]],
        pointers_to_wavefunctions[index[count_orbitals - 2]], size_order);
        t59_flag = true;
        }
    if (count_orbitals % 8 >= 1) {
        t60 = thread(&basis_set_calculations::Gradient_thread, this, pointers_to_Gradients[index[count_orbitals - 1]],
        pointers_to_wavefunctions[index[count_orbitals - 1]], size_order);
        t60_flag = true;
        }
    if (t54_flag == true) {
        t54.join();
        t54_flag = false;
        }
    if (t55_flag == true) {
        t55.join();
        t55_flag = false;
        }
    if (t56_flag == true) {
        t56.join();
        t56_flag = false;
        }
    if (t57_flag == true) {
        t57.join();
        t57_flag = false;
        }
    if (t58_flag == true) {
        t58.join();
        t58_flag = false;
        }
    if (t59_flag == true) {
        t59.join();
        t59_flag = false;
        }
    if (t60_flag == true) {
        t60.join();
        t60_flag = false;
        }
    // end of multithreading code
    // closed-shell basis set method optimalization code    
    for (i = 0; i < count_orbitals; i++)
        if (restriction == true and (spins[spin_paired[index[i]]] == -0.5)
        and (bonding[index[i]] == -1 or Z[index[i]] == Z[bonding[index[i]]])) // for restricted basis set method
            efective_radius_base_array[spin_paired[index[i]]] = efective_radius_base_array[index[i]];
    // end of closed-shell basis set method optimalization code
    if (alocate == true)
        {
        for (i = 0; i < count_electrons; i++) // copy values to vectors for restricted and unrestricted method
            {
            atom_wavefunctions->wavefunctions.push_back(pointers_to_wavefunctions[i]);
            atom_wavefunctions->probabilities.push_back(pointers_to_probabilities[i]);
            atom_wavefunctions->efective_radius_base.push_back(efective_radius_base_array[i]);
            atom_wavefunctions->Gradients.push_back(pointers_to_Gradients[i]);
            }
        }
    return(0);
    }
template <typename T>
T basis_set_calculations<T>::String_to_list_electrons(string UI_input, unsigned int size_order,
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
    string atoms_string[max_atoms];
    string coordinates_string[max_atoms];
    string bonds_string[max_atoms];
    string potentials_string[max_atoms];
    
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
        for (i = 0; i < count_atoms; i++) // sumarizing wavefunctions
            Sum_atomic_wavefunctions(&results, &wavefunctions[i]);
        }
    try {
        nuclear_atraction_integral_matrix = new T[matrix_order * matrix_order];
        nucleuses_atractions = new T[matrix_order * matrix_order];
        nucleuses_distances = new T[matrix_order * matrix_order];
        coulombic_integral_matrix = new T[matrix_order * matrix_order];
        overlap_integral_matrix = new T[matrix_order * matrix_order];
        overlap_efective_lenght_integral_matrix = new T[matrix_order * matrix_order];
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
        if (overlap_efective_lenght_integral_matrix != nullptr){
            delete[] overlap_efective_lenght_integral_matrix;
            overlap_efective_lenght_integral_matrix = nullptr;
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
T basis_set_calculations<T>::Calculate(unsigned int max_iterations, T minimal_fidelity, unsigned int size_order, bool alocate,
vector<T>* values, vector<T>* spin_density_vector,  vector<T>* spin_values)
    {
    unsigned int i, j;
    unsigned int matrix_order;
    T Hamiltonian, last_Hamiltonian;
    
    last_Hamiltonian = 0;
    matrix_order = results.n.size();
    memset(basis_set_matrix, 0, matrix_order * matrix_order);
    memset(corr_basis_set_matrix, 0, matrix_order * matrix_order);
    memset(spin_density_matrix, 0, matrix_order * matrix_order);
    if (alocate == true)
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
        Calculate_spin_density_matrix(overlap_integral_matrix, spin_density_matrix, matrix_order, &results);
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
            Calculate_resonance_integral_matrix(overlap_integral_matrix, overlap_efective_lenght_integral_matrix,
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
    Calculate_spin_density_matrix(overlap_integral_matrix, spin_density_matrix, matrix_order, &results);
    Solve_spin_density_matrix(spin_density_matrix, matrix_order, spin_values);
    for (i = 0; i < (matrix_order * matrix_order); i++)
        spin_density_vector->push_back(spin_density_matrix[i]);
    
    return(Hamiltonian);
    }
template <typename T>
T basis_set_calculations<T>::Clear()
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
    if (overlap_efective_lenght_integral_matrix != nullptr) {
        delete[] overlap_efective_lenght_integral_matrix;
        overlap_efective_lenght_integral_matrix = nullptr;
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
    results.wavefunction_coefficients.clear();
    results.wavefunction_lenght_multipliers.clear();
    results.efective_radius_base.clear();
    results.spins.clear();
    results.spin_paired.clear();
    results.bonding.clear();
    results.antibonding.clear();
    results.pi_bonding.clear();
    results.wavefunction_constraints.clear();
    results.n.clear();
    results.l.clear();
    results.m.clear();
    results.charge.clear();
    results.count_electrons.clear();
    results.reduced_Z.clear();
    results.Z.clear();
    results.electron_numbers.clear();
    results.x.clear();
    results.y.clear();
    results.z.clear();
    small_results.lenghts.clear();
    small_results.relative_lenghts.clear();
    small_results.wavefunctions.clear();
    small_results.probabilities.clear();
    small_results.electron_numbers.clear();
    small_results.lenght_orders.clear();
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
basis_set_calculations<T>::~basis_set_calculations(){
Clear();}
template class basis_set_calculations<double>;
/*
Author of this source code Ing. Pavel Florian Ph.D. licensed this source code under the the Apache License:
Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/
*/
