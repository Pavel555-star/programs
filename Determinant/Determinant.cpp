#include <iostream>
#include <string>
using namespace std;
template <typename T>
T determinant(int order, T* pointer) // Výpočet determinantu matice (spolehlivý pro všechny matice do 8 * 8 bytů čísel) při použití datového typu long
{
int i;
int j;
T det = 0;
T line = 1;

switch (order)
    {
    case 0:
        return(0); // Pro nesprávná zadání vrátí hodnotu 0
    case 1:
        return(pointer[0]); // Pro matice 1. řádu vrátí hodnotu prvku
    case 2:
        return((pointer[0]*pointer[3])-(pointer[1]*pointer[2])); // Pro matice 2. řádu vrátí determinant matice 2. řádu
    case 3:{ // Výpočet Sarusovým pravidlem
        T matrix[3][5];
        for (j = 0; j < 3; j++) // Překopírování matice do pomocné matice
            for (i = 0; i < 3; i++)
                matrix[i][j] = pointer[i+3*j];
        for (j = 0; j < 2; j++) // Naplnění pomocných řádků
            for (i = 0; i < 3; i++)
            matrix[i][j + 3] = matrix[i][j];
        for (j = 0; j < 3; j++) // Přičtení kladných řádků
            {
            for (i =0; i < 3; i++)
                line = line * matrix[i][i + j];
            det = det + line;
            line = 1;
                }      
        for (j = 0; j < 3; j++) // Odečtení záporných řádků
            {
            for (i =0; i < 3; i++)
                line = line * matrix[2 - i][i + j];
            det = det - line;
            line = 1;
                }
        return det;
    }
    }
if ((order >= 4) and (order <= 8)) // Výpočet rozvojem dle sloupce
{
    int x;
    int subi;
    int subj;
    int prefactor;
    T* submatrix = new T[(order-1) * (order-1)];
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
                    prefactor =  -1;
                det = det + (prefactor * pointer[order * x] * determinant(order-1, submatrix));
                }
    }
    delete submatrix;
    return det;
    }
    
if (order >= 9) // Výpočet Gaussovou eliminací do trojúhelníkového tvaru
{
    int radek;
    int k;
    int sign = 1;
    T den;
    T n;
    T d;
    T* buffer = new T[order * order];
    T* denominator = new T[order * order];
    T* temp1 = new T[order];
    T* temp2 = new T[order];
    for (i = 0; i < order; i++)
    {
        for (j = 0; j < order; j++)
        {
        buffer[i + (order * j)] = pointer[i + (order * j)]; // Překopírování matice k úpravám
        denominator[i + (order * j)] = 1; // Naplnění matice jmenovatelů jedničkami
        }
    }
    for (i = 0; i < order; i++)
    {
        for (j = i; j < order; j++)
        {
            if (buffer[i + (order * j)] != 0) // Nalezení prvního řádku s nenulovou hodnotou v daném sloupci
            {
            radek = j;
            break;
            }
            if ((j == (order-1)) and (buffer[i + (order * j)] == 0)) // Kontrola nuly
              return 0;
        }
        for (j = i; j < order; j++)
        {
            if ((radek != j) and (buffer[i + (order * j)] != 0)) // Vlastní přičítání nebo odečítání řádků
            {
                n = buffer[i + (order * j)];
                d = buffer[i + (order * radek)];
                if ((n % d) == 0)
                {
                n = n/d;
                d = 1;
                    for (k = 0; k < order; k++)
                        {
                        buffer[k + (order * j)] = buffer[k + (order * j)] - (buffer[k + (order * radek)] * n);
                        }
                }
                else
                    for (k = 0; k < order; k++)
                    {
                    buffer[k + (order * j)] = (d * buffer[k + (order * j)]) - (buffer[k + (order * radek)] * n);
                    denominator[k + (order * j)] = denominator[k + (order * j)] * d;
                    }
            }
            if (i != radek) // Prohození řádků
                {
                sign = sign * -1;
                for (k = 0; k < order; k++)
                    {
                    temp1[k] = buffer[k + (order * radek)];
                    temp2[k] = denominator[k + (order * radek)];
                    buffer[k + (order * radek)] = buffer[k + (order * i)];
                    denominator[k + (order * radek)] = denominator[k + (order * i)];
                    buffer[k + (order * i)] = temp1[k];
                    denominator[k + (order * i)] = temp2[k];
                    }
                }
        }
    }
    det = 1;
    den = 1;
    for (i = 0; i < order; i++)
    {
        det = det * buffer[i + (order * i)];
        den = den * denominator[i + (order * i)];
        if ((det % den) == 0)
        {
        det = det/den;
        den =1;
        }
    }
    det = (det/den)*sign;

delete denominator;
delete buffer;
delete temp1;
delete temp2;
return det;
}
}
int main(void) {

int order, j, i;
long* matrix;

cout << "Zadejte řád matice / Enter the order of matrix:\n";
cin >> order;
matrix = new long[(order) * (order)];

cout << "Zadejte elementy matice / Enter the elements of matrix:\n";
for (j = 0; j < order; j++)
   for (i = 0; i < order; i++)
      cin >> matrix[i + (order * j)];

cout << "Zadali jste matici / You are entered the matrix: \n";     
   for (j = 0; j < order; j++) {  
      for (i = 0; i < order; i++)
         cout << matrix[i + (order * j)] <<" ";
      cout<<endl;
      }
cout<<"Determinant matice je / Determinant of the matrix is: "<< determinant(order, matrix) << " ";

return 0;} // Autor Ing. Pavel Florián Ph.D. dovoluje šíření tohoto programu jen s uvedením svého jména.
