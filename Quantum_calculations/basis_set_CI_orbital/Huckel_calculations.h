#include <iostream>
#include <string>
#include <cmath>
#include <thread> // compiler parameters: -pthread
#include <vector>
#include <bits/stdc++.h>
using namespace std;

#ifndef HUCKEL_CALCULATIONS_H
#define HUCKEL_CALCULATIONS_H

template <typename T>
class Huckel_calculations
{
public:
    vector<T> determinants;
    vector<T> spectra;
    vector<T> coefficients;
    vector<T> sequence;
    
    bool cyclic;
    unsigned int determinant_exception_handle = 0;
    
    inline T Determinant(unsigned int order, T* matrix_pointer, T* buffer, T* denominator, T* temp1, T* temp2);
    T Huckel_Determinant_check(unsigned int order, T* matrix);
    T Huckel_Determinant_set(unsigned int order, T* pointer,unsigned int count, T min, T step, T* output_values);
    T Huckel_Determinant_solver(unsigned int order, T* matrix);
    T Huckel_Spectra_solver(); // Solving of energy differences between occupied and unoccupied energy levels
    T Huckel_Wavefunction_Coefficient_solver(unsigned int order, T* matrix);
};
template <typename T>
inline T Huckel_calculations<T>::Determinant(unsigned int order, T* pointer, T* buffer, T* denominator, T* temp1, T* temp2) 
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
T Huckel_calculations<T>::Huckel_Determinant_check(unsigned int order, T* matrix)
    {
    unsigned int i;
    unsigned int j;
    unsigned int count;
    
    if (order < 2)
        return(-1);
    
    for (i = 0; i < order; i++) // Check for number of non diagonal members > 0 between 1-3
    {
        count =0;
        for (j = 0; j < order; j++)
        {
            if ((matrix[(i * order) + j] > 0) and (i != j))
            {
                count++;
                if (matrix[(i * order) + j] == (matrix[(j * order) + i] == 0))
                    {
                    return(-1);
                    }
            }
        }
        if ((count < 1) or (count > 3))
        {
        return(-1);
        }        
    }
    return(0);
    }
template <typename T>    
T Huckel_calculations<T>::Huckel_Determinant_set(unsigned int order, T* pointer,unsigned int count,
T min, T step, T* output_values)
    {
    unsigned int i;
    unsigned int j;
    
    if (order < 2)
        return(-1);
    
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
        if (order > 32)
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
    {
    auxiliar_matrix[i] = pointer[i];
    }
    for (i = 0; i < order; i++) // Copying the diagonal of the matrix
    {
    diagonal[i] = auxiliar_matrix[i * (order + 1)];
    }   
    for (i = 0; i < count; i++) // Calculating the determinant of te values array
    {
        for (j = 0; j < order; j++)
            {
            auxiliar_matrix[j * (order + 1)] =  value + diagonal[j];
            }
        output_values[i] = Determinant(order, auxiliar_matrix, buffer, denominator, temp1, temp2);
    value = value + step;
    }
    delete[] diagonal;
    delete[] auxiliar_matrix;
    
    if (order > 32)
        {
        delete[] buffer;
        delete[] denominator;
        delete[] temp1;
        delete[] temp2;
        }
    
    return(0);
    }
template <typename T>      
T Huckel_calculations<T>::Huckel_Determinant_solver(unsigned int order, T* matrix)
    {
    unsigned int i, j;
    unsigned int r_m;
    
    T* determinant_values = nullptr;
    bool extended = false;
    T determinant_max = 0;
    T row = 0;
    T row_max = 0;
    unsigned int range_multiplier;
    
    try {
        determinant_values = new T[64000];
        } // Calculating of determinant values array
    catch (int){
        determinant_exception_handle = 1;
        return (-1);
        }
    
    
    for (i = 0; i < order; i++)
        {
        row = 0;
        for (j = 0; j < order; j++)
            {
            row = row + matrix[(i * order) + j];
            }
        if (row > row_max)
            row_max = row;
        }
    if (row_max < 3)
        r_m = 1;
    else
        if (row_max < 6)
            r_m = 2;
        else
            r_m = 8;
    
    thread t1(&Huckel_calculations<T>::Huckel_Determinant_set, this, order, matrix, 8000, -4.00 * r_m, 0.000125 * r_m,
    determinant_values);
    thread t2(&Huckel_calculations<T>::Huckel_Determinant_set, this, order, matrix, 8000, -3.00 * r_m, 0.000125 * r_m,
    determinant_values +  8000);
    thread t3(&Huckel_calculations<T>::Huckel_Determinant_set, this, order, matrix, 8000, -2.00 * r_m, 0.000125 * r_m,
    determinant_values + 16000);
    thread t4(&Huckel_calculations<T>::Huckel_Determinant_set, this, order, matrix, 8000, -1.00 * r_m, 0.000125 * r_m,
    determinant_values + 24000);
    thread t5(&Huckel_calculations<T>::Huckel_Determinant_set, this, order, matrix, 8000,  0.00 * r_m, 0.000125 * r_m,
    determinant_values + 32000);
    thread t6(&Huckel_calculations<T>::Huckel_Determinant_set, this, order, matrix, 8000,  1.00 * r_m, 0.000125 * r_m,
    determinant_values + 40000);
    thread t7(&Huckel_calculations<T>::Huckel_Determinant_set, this, order, matrix, 8000,  2.00 * r_m, 0.000125 * r_m,
    determinant_values + 48000);
    thread t8(&Huckel_calculations<T>::Huckel_Determinant_set, this, order, matrix, 8000,  3.00 * r_m, 0.000125 * r_m,
    determinant_values + 56000);
    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();
    t7.join();
    t8.join();
      
    determinants.clear();
    for (i = 1; i < 63999; i++)
        if (determinant_values[i] * determinant_values[i] > determinant_max)
            determinant_max = determinant_values[i] * determinant_values[i];
            
    for (i = 1; i < 63999; i++)
        {
        if ((determinant_values[i] * determinant_values[i] <= determinant_max * 0.1)
        and (determinant_values[i] * determinant_values[i] <  determinant_values[i + 1] * determinant_values[i + 1])
        and (determinant_values[i] * determinant_values[i] < determinant_values[i - 1] * determinant_values[i - 1]))
            determinants.push_back((T(i)/(8000/r_m)) -(4.00 * r_m));
        }
    delete(determinant_values);
    return(0);
    }
    
template <typename T>    
T Huckel_calculations<T>::Huckel_Spectra_solver() 
    { // Solving of energy differences between occupied and unoccupied energy levels
    unsigned int i;
    unsigned int j;
    unsigned int size;
    
    T S; // Overlap integral
    T* spectra_matrix;
        
    spectra.clear();
    size = determinants.size();
    if (S > 1)
        S = 1;
    if (S < 0)
        S = 0;
    if (size > 2)
        {
        spectra_matrix = new T[size * size];
        
        for (i = 0; i < size * size; i++) // Filling matrix for spectra by zeros
            spectra_matrix[i] = 0;
        
        for (i = 0; i < (size/2 + size % 2); i++) // Calculating energies of transitions
            for (j = (size/2); j < size; j++)
                spectra_matrix[i + j * (size/2 + size % 2)] = determinants[j] - determinants[i];
            
            for (i = 0; i < size * size; i++)
                if (spectra_matrix[i] > 0)
                    spectra.push_back(spectra_matrix[i]);
                
            sort(spectra.begin(), spectra.end()); // Sorting values and copying into spectra vector
        
        delete[] spectra_matrix;
        }
    else
        return(-1);
    return(0);
    }
template <typename T>    
T Huckel_calculations<T>::Huckel_Wavefunction_Coefficient_solver(unsigned int order, T* matrix)
    {
    unsigned int i, j, k, l;
    unsigned int cycle_max;
    T auxiliar_matrix_max;
    auxiliar_matrix_max = 0;
    if (order < 2)
        return(-1);
        
    unique_ptr<unsigned int[]> bond_count(new unsigned int[order]);
    unique_ptr<unsigned int[]> path(new unsigned int[order]);
    unique_ptr<unsigned int[]> non_path(new unsigned int[order]);
    unique_ptr<bool[]> processed_lines(new bool[order]);
    unique_ptr<T[]> auxiliar_matrix(new T[order * order]);
    unique_ptr<T[]> auxiliar_matrix_2(new T[order * order]);
    unique_ptr<T[]> coefficient_matrix(new T[order * order]);
    unique_ptr<T[]> circle_positions(new T[order]);
    
    unsigned int path_size = 0;
    unsigned int bond_counts = 0;
    unsigned int erased_count = 1;
    unsigned int erased_counts = 0;
    int first_one_bonding = -1;
    int first_two_bonding = -1;
    int count_item_to_process;
    
    T norm_coefficient = 0;
    cycle_max = order * order;
    coefficients.clear(); // Clear coefficient values vector  
              
    for (i = 0; i < (order * order); i++) 
        {
        auxiliar_matrix[i] = matrix[i]; // Copying of matrix
        auxiliar_matrix_2[i] = matrix[i];
        coefficient_matrix[i] = 0; // Filling coefficient matrix by 0
        }
    for (i = 0; i < order; i++) // Filling the diagonal of auxilliary matrix, bond_count path and circle_positions array by zeros
        {
        auxiliar_matrix[i * (order + 1)] = 0;
        bond_count[i] = 0;
        path[i] = 0;
        circle_positions[i] = 0;
        }
    for (i = 0; i < (order * order); i++) // search maximal value of auxiliar matrix cell
        if (auxiliar_matrix[i] > auxiliar_matrix_max)
            auxiliar_matrix_max = auxiliar_matrix[i];
    for (i = 0; i < (order * order); i++) // preparation of auxiliar matrix
        if (auxiliar_matrix[i] < auxiliar_matrix_max * 0.20 and i % order != i/order)
            auxiliar_matrix[i] = 0;
    
    for (i = 0; i < order; i++) // Detection of cyclicity 
        {                           // and searching for first one and two-bonding C atoms
        for (j = 0; j < order; j++)
            if (auxiliar_matrix[(i * order) + j] > 0)
                bond_count[i]++; 
                
        if ((bond_count[i] == 1) and (first_one_bonding == -1))
                first_one_bonding = i;
                
        if ((bond_count[i] == 2) and (first_two_bonding == -1))
                first_two_bonding = i;
                
        bond_counts += bond_count[i];           
        }
    auto Cyclic = [](unsigned int bond_count, unsigned int order) {return((bond_count/2) >= order); }; // Lambda function
    cyclic = Cyclic(bond_counts, order);
        
    if (cyclic == false and first_one_bonding != -1) // Path finding for non-cyclic hydrocarbons
        {
        path[0] = first_one_bonding;
        path_size = 1;
        for (i = 1; i < order; i++)
            {
            for (j = 0; j < order; j++)
                {
                if (i >= 2) // continue with new value of j with previous item
                    if (j == path[i - 2])
                        continue;
                        
                if ((auxiliar_matrix[(path[i-1] * order) + j] > 0)) // Choose next atom in string
                    break;
                }
            path[i] = j;
            path_size += 1;
            }
        }
    if (cyclic == true and first_two_bonding != -1) // Path finding for cyclic hydrocarbons
        {
        if (first_one_bonding != -1) // Delete linear chains
            {
            while (erased_count != 0 and cycle_max != 0)
                {
                erased_count = 0;
                for (i = 1; i < order; i++)
                    if (bond_count[i] == 1)
                        for (j = 0; j < order; j++)
                            if (auxiliar_matrix[(i * order) + j] > 0)
                                {
                                auxiliar_matrix[(i * order) + j] = 0;
                                auxiliar_matrix[(j * order) + i] = 0;
                                bond_count[i]--;
                                bond_count[j]--;
                                erased_count++;
                                erased_counts++;
                                }
                cycle_max--;
                }
            }
        path[0] = first_two_bonding; // Finding the path
        path_size = 1;
        for (i = 0; i < (order - erased_counts); i++)
            {
            for (j = 0; j < order; j++)
                {
                if (i > 0)
                    {
                    if ((auxiliar_matrix[(path[i] * order) + j] > 0) and ((j != path[i-1]) or i == 0) and (j != path[0]))
                        {  // Conditions of path selection
                        if (bond_count[j] < 3)
                            {
                            path[i + 1] = j;
                            path_size += 1;
                            break; 
                            }
                        if ((bond_count[j] == 3) and (path[i] == 2))
                            {
                            path[i + 1] = j;
                            path_size += 1;
                            break;
                            }
                        }
                    }
                else
                    {
                    if ((auxiliar_matrix[(path[i] * order) + j] > 0) and (j != path[0]))
                        {  // Conditions of path selection
                        if (bond_count[j] < 3)
                            {
                            path[i + 1] = j;
                            path_size += 1;
                            break; 
                            }
                        if ((bond_count[j] == 3) and (path[i] == 2))
                            {
                            path[i + 1] = j;
                            path_size += 1;
                            break;
                            }
                        }
                    }
                }
            }
        }
    for (i =0; i < path_size; i++)
        sequence.push_back(path[i]);
       
    for (i =0; i < order; i++) // Calculation of wavefunction coefficients
        {
        for (j = 0; j < order; j++)
            processed_lines[j] = false;
                
        count_item_to_process = order;
        for (j =0; j < order; j++) // Calcuation of coefficients for atoms in path
            {
            if ((cyclic == false) and (j < path_size)) //  linear hydrocarbons
                {
                coefficient_matrix[path[j] + (i * order)] = sin(M_PI * (i + 1) * (j + 1)/(order + 1));
                processed_lines[path[j]] = true;
                count_item_to_process--;
                }
            if ((cyclic == true)  and (j < path_size)) // cyclic hydrocarbons - atoms in path
                {
                if ((i % 2 == 0) or (i + 1 == order))
                    {
                    circle_positions[path[j]] = 2 * M_PI * ((i + 1)/2) * (j + 1)/(order);
                    coefficient_matrix[path[j] + (i * order)] = cos(circle_positions[j]);
                    }
                else
                    {
                    circle_positions[path[j]] = 2 * M_PI * ((i + 1)/2) * (j + 1)/(order);
                    coefficient_matrix[path[j] + (i * order)] = sin(circle_positions[j]);
                    }
                count_item_to_process--;
                processed_lines[path[j]] = true;
                }
            }
        if ((cyclic == true) and (order >= path_size))//cyclic hydrocarbons - search atoms out of path and its circle_positions 
            {
            l = 0;
            cycle_max = order * order;
            while (count_item_to_process > 0 and cycle_max != 0)
                {
                for (j = 0; j < order; j++)
                    if (processed_lines[j] == false)
                        {
                        if (l > (order - path_size)) // Filling of list of non-path atoms
                            {
                            non_path[l] = j;
                            l++;
                            }
                        for (k = 0; k < order; k++)
                            {
                            if ((auxiliar_matrix_2[k + (j * order)] > 0) and (processed_lines[k] == true))
                                {
                                circle_positions[j] = circle_positions[k] + (2 * M_PI * ((i + 1)/2)/double(order));
                                if (processed_lines[j] == false)
                                    {
                                    processed_lines[j] = true;
                                    non_path[l] = j;
                                    l++;
                                    }
                                }
                            }
                        if (processed_lines[j] == true)
                            {
                            count_item_to_process--;
                            if ((i % 2 == 0) or (i + 1 == order))
                                coefficient_matrix[j + (i * order)] = cos(circle_positions[j]);
                            else
                                coefficient_matrix[j + (i * order)] = sin(circle_positions[j]);
                            }
                        }
                cycle_max--;
                }
            } 
        for (j =0; j < order; j++) // Calculation of normalization coefficient
            norm_coefficient +=  coefficient_matrix[j + (i * order)] * coefficient_matrix[j + (i * order)];
                
        norm_coefficient = sqrt(norm_coefficient);
        for (j =0; j < order; j++) // Normalization of wavefunctions coefficients
            {
            coefficient_matrix[j + (i * order)] /= norm_coefficient;
            coefficients.push_back(coefficient_matrix[j + (i * order)]);
            }
        norm_coefficient = 0;
        }
    return(0);
    }
# endif // HUCKEL_CALCULATIONS_H
/*
Author of this source code Ing. Pavel Florian Ph.D. licensed this source code under the the Apache License:

Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/
*/
