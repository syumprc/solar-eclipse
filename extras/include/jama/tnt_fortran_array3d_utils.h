/*
*
* Template Numerical Toolkit (TNT)
*
*
*/


#ifndef TNT_FORTRAN_ARRAY3D_UTILS_H
#define TNT_FORTRAN_ARRAY3D_UTILS_H

#include <cstdlib>
#include <cassert>

namespace TNT
{


template <class T>
std::ostream& operator<<(std::ostream &s, const Fortran_Array3D<T> &A)
{
    int M=A.dim1();
    int N=A.dim2();
    int K=A.dim3();

    s << M << " " << N << " " << K << "\n";

    for (int i=1; i<=M; i++)
    {
        for (int j=1; j<=N; j++)
        {
			for (int k=1; k<=K; k++)
            	s << A(i,j,k) << " ";
			s << "\n";
        }
        s << "\n";
    }


    return s;
}

template <class T>
std::istream& operator>>(std::istream &s, Fortran_Array3D<T> &A)
{

    int M, N, K;

    s >> M >> N >> K;

	Fortran_Array3D<T> B(M,N,K);

    for (int i=1; i<=M; i++)
        for (int j=1; j<=N; j++)
			for (int k=1; k<=K; k++)
            	s >>  B(i,j,k);

	A = B;
    return s;
}


template <class T>
Fortran_Array3D<T> operator+(const Fortran_Array3D<T> &A, const Fortran_Array3D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();
	int p = A.dim3();

	if (B.dim1() != m ||  B.dim2() != n || B.dim3() != p )
		return Fortran_Array3D<T>();

	else
	{
		Fortran_Array3D<T> C(m,n,p);

		for (int i=1; i<=m; i++)
			for (int j=1; j<=n; j++)
				for (int k=1; k<=p; k++)
				C(i,j,k) = A(i,j,k)+ B(i,j,k);

		return C;
	}
}


template <class T>
Fortran_Array3D<T> operator-(const Fortran_Array3D<T> &A, const Fortran_Array3D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();
	int p = A.dim3();

	if (B.dim1() != m ||  B.dim2() != n || B.dim3() != p )
		return Fortran_Array3D<T>();

	else
	{
		Fortran_Array3D<T> C(m,n,p);

		for (int i=1; i<=m; i++)
			for (int j=1; j<=n; j++)
				for (int k=1; k<=p; k++)
				C(i,j,k) = A(i,j,k)- B(i,j,k);

		return C;
	}
}


template <class T>
Fortran_Array3D<T> operator*(const Fortran_Array3D<T> &A, const Fortran_Array3D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();
	int p = A.dim3();

	if (B.dim1() != m ||  B.dim2() != n || B.dim3() != p )
		return Fortran_Array3D<T>();

	else
	{
		Fortran_Array3D<T> C(m,n,p);

		for (int i=1; i<=m; i++)
			for (int j=1; j<=n; j++)
				for (int k=1; k<=p; k++)
				C(i,j,k) = A(i,j,k)* B(i,j,k);

		return C;
	}
}


template <class T>
Fortran_Array3D<T> operator/(const Fortran_Array3D<T> &A, const Fortran_Array3D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();
	int p = A.dim3();

	if (B.dim1() != m ||  B.dim2() != n || B.dim3() != p )
		return Fortran_Array3D<T>();

	else
	{
		Fortran_Array3D<T> C(m,n,p);

		for (int i=1; i<=m; i++)
			for (int j=1; j<=n; j++)
				for (int k=1; k<=p; k++)
				C(i,j,k) = A(i,j,k)/ B(i,j,k);

		return C;
	}
}


template <class T>
Fortran_Array3D<T>& operator+=(Fortran_Array3D<T> &A, const Fortran_Array3D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();
	int p = A.dim3();

	if (B.dim1() == m &&  B.dim2() == n && B.dim3() == p )
	{
		for (int i=1; i<=m; i++)
			for (int j=1; j<=n; j++)
				for (int k=1; k<=p; k++)
					A(i,j,k) += B(i,j,k);
	}

	return A;
}


template <class T>
Fortran_Array3D<T>& operator-=(Fortran_Array3D<T> &A, const Fortran_Array3D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();
	int p = A.dim3();

	if (B.dim1() == m &&  B.dim2() == n && B.dim3() == p )
	{
		for (int i=1; i<=m; i++)
			for (int j=1; j<=n; j++)
				for (int k=1; k<=p; k++)
					A(i,j,k) -= B(i,j,k);
	}

	return A;
}


template <class T>
Fortran_Array3D<T>& operator*=(Fortran_Array3D<T> &A, const Fortran_Array3D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();
	int p = A.dim3();

	if (B.dim1() == m &&  B.dim2() == n && B.dim3() == p )
	{
		for (int i=1; i<=m; i++)
			for (int j=1; j<=n; j++)
				for (int k=1; k<=p; k++)
					A(i,j,k) *= B(i,j,k);
	}

	return A;
}


template <class T>
Fortran_Array3D<T>& operator/=(Fortran_Array3D<T> &A, const Fortran_Array3D<T> &B)
{
	int m = A.dim1();
	int n = A.dim2();
	int p = A.dim3();

	if (B.dim1() == m &&  B.dim2() == n && B.dim3() == p )
	{
		for (int i=1; i<=m; i++)
			for (int j=1; j<=n; j++)
				for (int k=1; k<=p; k++)
					A(i,j,k) /= B(i,j,k);
	}

	return A;
}


} // namespace TNT

#endif
