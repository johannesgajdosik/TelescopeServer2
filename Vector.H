/*
 * Author and Copyright
 * Johannes Gajdosik, 2008
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <iostream>

#include <cmath>

template <typename T,int dim>
class Vector {
public:
    // type introspection
  typedef T Scal;
  enum {dimension = dim};

    // construction
  Vector(void) {}
  explicit Vector(const T &x0)
    {c[0]=x0;}
  Vector(const T &x0,const T &x1)
    {c[0]=x0;c[1]=x1;}
  Vector(const T &x0,const T &x1,const T &x2)
    {c[0]=x0;c[1]=x1;c[2]=x2;}
  Vector(const T &x0,const T &x1,const T &x2,const T &x3)
    {c[0]=x0;c[1]=x1;c[2]=x2;c[3]=x3;}
  Vector(const T &x0,const T &x1,const T &x2,const T &x3,const T &x4)
    {c[0]=x0;c[1]=x1;c[2]=x2;c[3]=x3;c[4]=x4;}
  Vector(const T &x0,const T &x1,const T &x2,const T &x3,const T &x4,const T &x5)
    {c[0]=x0;c[1]=x1;c[2]=x2;c[3]=x3;c[4]=x4;c[5]=x5;}

  Vector(const T v[])
    {for (int i=0;i<dim;i++) c[i] = v[i];}
  const Vector &operator=(const T v[])
    {for (int i=0;i<dim;i++) {c[i] = v[i];} return *this;}

  template <typename U> Vector(const Vector<U,dim> &v)
    {for (int i=0;i<dim;i++) c[i] = v[i];}
  template <typename U> const Vector &operator=(const Vector<U,dim> &v)
    {for (int i=0;i<dim;i++) {c[i] = v[i];} return *this;}

    // member access
  const T &operator[](int i) const {return c[i];}
  T &operator[](int i) {return c[i];}
  operator const T*(void) const {return c;}
  operator T*(void) {return c;}
  void setZero(void) {for (int i=0;i<dim;i++) c[i] = 0;}
  bool nonZero(void) const
    {for (int i=0;i<dim;i++) {if (c[i]) return true;} return false;}

    // arithmetic
  template <typename U> const Vector &operator+=(const Vector<U,dim> &v)
    {for (int i=0;i<dim;i++) {c[i] += v[i];} return *this;}
  template <typename U> const Vector &operator-=(const Vector<U,dim> &v)
    {for (int i=0;i<dim;i++) {c[i] -= v[i];} return *this;}
  Vector operator-(void) const
    {Vector v;for (int i=0;i<dim;i++) {v[i] = -c[i];} return v;}
  const Vector &negate(void)
    {for (int i=0;i<dim;i++) {c[i] = -c[i];} return *this;}

  const Vector &operator*=(const T &x)
    {for (int i=0;i<dim;i++) {c[i] *= x;} return *this;}

  template <typename U> const Vector &operator*=(const U &x)
    {for (int i=0;i<dim;i++) {c[i] *= x;} return *this;}
  template <typename U> const Vector &operator/=(const U &x)
    {for (int i=0;i<dim;i++) {c[i] /= x;} return *this;}

    // geometric functions
  inline T length2(void) const;
protected:
  T c[dim];
};

template <typename T,typename U>struct Result {};
template <> struct Result<double,double> {typedef double Type;};
template <> struct Result<double,float> {typedef double Type;};
template <> struct Result<float,double> {typedef double Type;};
template <> struct Result<double,int> {typedef double Type;};
template <> struct Result<int,double> {typedef double Type;};
template <> struct Result<float,float> {typedef float Type;};
template <> struct Result<float,int> {typedef float Type;};
template <> struct Result<int,float> {typedef float Type;};

template <typename T,typename U>struct MulResult {};
template <> struct MulResult<double,double> {typedef double Type;};
template <> struct MulResult<double,float> {typedef double Type;};
template <> struct MulResult<float,double> {typedef double Type;};
template <> struct MulResult<double,int> {typedef double Type;};
template <> struct MulResult<int,double> {typedef double Type;};
template <> struct MulResult<float,float> {typedef float Type;};
template <> struct MulResult<float,int> {typedef float Type;};
template <> struct MulResult<int,float> {typedef float Type;};

class GmpFloat;
template <> struct Result<GmpFloat,GmpFloat> {typedef GmpFloat Type;};
template <> struct MulResult<GmpFloat,GmpFloat> {typedef GmpFloat Type;};

template <int dim> struct MulResult<Vector<double,dim>,double> {
  typedef Vector<double,dim> Type;
};

template <int dim> struct MulResult<double,Vector<double,dim> > {
  typedef Vector<double,dim> Type;
};

// how to do this?
//template <>
//template <typename T,int dim,typename U>
//struct MulResult<Vector<T,dim>,U> {
//  typedef Vector<
//            MulResult<
//              T,
//              U
//            >::Type,
//            dim
//          > Type;
//};


  // linearity
template <typename T,typename U,int dim> inline
Vector<typename Result<T,U>::Type,dim> operator+(const Vector<T,dim> &v,
                                                 const Vector<U,dim> &w) {
  return Vector<typename Result<T,U>::Type,dim>(v) += w;
}

template <typename T,typename U,int dim> inline
Vector<typename Result<T,U>::Type,dim> operator-(const Vector<T,dim> &v,
                                                 const Vector<U,dim> &w) {
  return Vector<typename Result<T,U>::Type,dim>(v) -= w;
}

template <typename T,typename U,int dim> inline
Vector<typename MulResult<T,U>::Type,dim> operator*(const T &x,
                                                    const Vector<U,dim> &v) {
  return Vector<typename MulResult<T,U>::Type,dim>(v) *= x;
}

template <typename T,typename U,int dim> inline
Vector<typename MulResult<T,U>::Type,dim> operator*(const Vector<T,dim> &v,
                                                    const U &x) {
  return Vector<typename MulResult<T,U>::Type,dim>(v) *= x;
}

template <typename T,typename U,int dim> inline
Vector<typename MulResult<T,U>::Type,dim> operator/(const Vector<T,dim> &v,
                                                    const U &x) {
  return Vector<typename MulResult<T,U>::Type,dim>(v) /= x;
}



  // inner product
template <typename T,typename U,int dim> inline
typename Result<T,U>::Type operator|(const Vector<T,dim> &v,
                                     const Vector<U,dim> &w) {
  typename Result<T,U>::Type x(0);
  for (int i=0;i<dim;i++) x += v[i]*w[i];
  return x;
}

template <typename T,int dim>
T Vector<T,dim>::length2(void) const {
  return operator|(*this,*this);
}



template <typename T,typename U> inline
Vector<typename Result<T,U>::Type,3> operator&(const Vector<T,3> &v,
                                               const Vector<U,3> &w) {
  Vector<typename Result<T,U>::Type,3> x;
  x[0] = v[1]*w[2] - v[2]*w[1];
  x[1] = v[2]*w[0] - v[0]*w[2];
  x[2] = v[0]*w[1] - v[1]*w[0];
  return x;
}

template <typename T>
Vector<T,3> PolarToRect(T ra,T dec) {
  const T cd = cos(dec);
  return Vector<T,3>(cos(ra)*cd,sin(ra)*cd,sin(dec));
}

template <typename T>
void RectToPolar(const Vector<T,3> &x,T &ra,T &dec) {
  ra = atan2(x[1],x[0]); // well defined even if both are 0
  if (ra < 0.0) ra += 2.0*M_PI;
  dec = atan2(x[2],sqrt(x[0]*x[0]+x[1]*x[1]));
}


  // IO
template <typename T,int dim> inline
std::ostream &operator<<(std::ostream &os,const Vector<T,dim> &v) {
  os << '[';
  if (dim > 0) {
    os << v[0];
    for (int i=1;i<dim;i++) os << ',' << v[i];
  }
  return os << ']';
}

template <typename T,int dim> inline
std::istream &operator>>(std::istream &is,Vector<T,dim> &v) {
  char ch;
  if (!(is >> ch)) return is;
  if (ch != '[') {
    is.setstate(std::ios::failbit|std::ios::badbit);
    return is;
  }
  if (dim > 0) {
    if (!(is >> v[0])) return is;
    for (int i=1;i<dim;i++) {
      if (!(is >> ch)) return is;
      if (ch != ',') {
        is.setstate(std::ios::failbit|std::ios::badbit);
        return is;
      }
      if (!(is >> v[i])) return is;
    }
  }
  if (!(is >> ch)) return is;
  if (ch != ']') is.setstate(std::ios::failbit|std::ios::badbit);
  return is;
}









template <typename T,int m,int n>
class Matrix {
public:
    // type introspection
  typedef T Scal;
  enum {cols = m,rows = n,dim = m*n};

    // construction
  Matrix(void) {}
  explicit Matrix(const T &x00)
    {x[0]=x00;}
  Matrix(const T &x00,const T &x01,const T &x10,const T &x11)
    {x[0]=x00;x[1]=x01;x[2]=x10;x[3]=x11;}
  Matrix(const T &x0,const T &x1,const T &x2,
         const T &x3,const T &x4,const T &x5,
         const T &x6,const T &x7,const T &x8) {
    x[0]=x0;x[1]=x1;x[2]=x2;
    x[3]=x3;x[4]=x4;x[5]=x5;
    x[6]=x6;x[7]=x7;x[8]=x8;
  }

  Matrix(const T a[])
    {for (int i=0;i<dim;i++) x[i] = a[i];}
  const Matrix &operator=(const T a[])
    {for (int i=0;i<dim;i++) {x[i] = a[i];} return *this;}

  template <typename U> Matrix(const Matrix<U,m,n> &a)
    {for (int i=0;i<dim;i++) x[i] = a[i];}
  template <typename U> const Matrix &operator=(const Matrix<U,m,n> &a)
    {for (int i=0;i<dim;i++) {x[i] = a[i];} return *this;}
  template <typename U> const Matrix &operator=(const U *a)
    {for (int i=0;i<dim;i++) {x[i] = a[i];} return *this;}

    // member access
    // TODO: rows/columns are mixed up. see
    // http://en.wikipedia.org/wiki/Index_notation
  const T &operator()(int i,int j) const {return x[i+m*j];}
  T &operator()(int i,int j) {return x[i+m*j];}
  Vector<T,rows> col(int i) const {
    Vector<T,rows> rval;
    for (int j=0;j<rows;j++) rval[j] = x[i+m*j];
    return rval;
  }
  Vector<T,cols> row(int j) const {
    return Vector<T,cols>(x+m*j);
  }
  operator const T*(void) const {return x;}
  operator T*(void) {return x;}
  void setZero(void) {for (int i=0;i<dim;i++) x[i] = 0;}
  bool nonZero(void) const
    {for (int i=0;i<dim;i++) {if (x[i]) return true;} return false;}

    // arithmetic
  template <typename U> const Matrix &operator+=(const Matrix<U,m,n> &a)
    {for (int i=0;i<dim;i++) {x[i] += a.x[i];} return *this;}
  template <typename U> const Matrix &operator-=(const Matrix<U,m,n> &a)
    {for (int i=0;i<dim;i++) {x[i] -= a.x[i];} return *this;}
  Matrix operator-(void) const
    {Matrix a;for (int i=0;i<dim;i++) {a[i] = -x[i];} return a;}
  template <typename U> const Matrix &operator*=(const U &u)
    {for (int i=0;i<dim;i++) {x[i] *= u;} return *this;}
  template <typename U> const Matrix &operator/=(const U &u)
    {for (int i=0;i<dim;i++) {x[i] /= u;} return *this;}

  Matrix<T,n,m> transpose(void) const {
    Matrix<T,n,m> a;
    for (int i=0;i<m;i++) for (int j=0;j<n;j++) a(j,i) = operator()(i,j);
    return a;
  }
protected:
  T x[dim];
};

  // linearity
template <typename T,typename U,int m,int n> inline
Matrix<typename Result<T,U>::Type,m,n> operator+(const Matrix<T,m,n> &a,
                                                 const Matrix<U,m,n> &b) {
  return Matrix<typename Result<T,U>::Type,m,n>(a) += b;
}

template <typename T,typename U,int m,int n> inline
Matrix<typename Result<T,U>::Type,m,n> operator-(const Matrix<T,m,n> &a,
                                                 const Matrix<U,m,n> &b) {
  return Matrix<typename Result<T,U>::Type,m,n>(a) -= b;
}

template <typename T,typename U,int m,int n> inline
Matrix<typename MulResult<T,U>::Type,m,n> operator*(const T &x,
                                                    const Matrix<U,m,n> &a) {
  return Matrix<typename MulResult<T,U>::Type,m,n>(a) *= x;
}

template <typename T,typename U,int m,int n> inline
Matrix<typename MulResult<T,U>::Type,m,n> operator*(const Matrix<T,m,n> &a,
                                                    const U &x) {
  return Matrix<typename MulResult<T,U>::Type,m,n>(a) *= x;
}

template <typename T,typename U,int m,int n> inline
Matrix<typename MulResult<T,U>::Type,m,n> operator/(const Matrix<T,m,n> &a,
                                                    const U &x) {
  return Matrix<typename MulResult<T,U>::Type,m,n>(a) /= x;
}

  // matrix/vector product
template <typename T,typename U,int l,int m,int n> inline
Matrix<typename Result<T,U>::Type,n,l> operator*(const Matrix<T,m,l> &a,
                                                 const Matrix<U,n,m> &b) {
  Matrix<typename Result<T,U>::Type,n,l> c;
  for (int i=0;i<l;i++) for (int k=0;k<n;k++) {
    c(k,i) = a(0,i)*b(k,0);
    for (int j=1;j<m;j++) c(k,i) += a(j,i)*b(k,j);
  }
  return c;
}

template <typename T,typename U,int l,int m> inline
Vector<typename Result<T,U>::Type,l> operator*(const Matrix<T,m,l> &a,
                                               const Vector<U,m> &b) {
  Vector<typename Result<T,U>::Type,l> c;
  for (int i=0;i<l;i++) {
    c[i] = a(0,i)*b[0];
    for (int j=1;j<m;j++) c[i] += a(j,i)*b[j];
  }
  return c;
}

template <typename T,typename U,int m,int n> inline
Vector<typename Result<T,U>::Type,n> operator*(const Vector<T,m> &a,
                                               const Matrix<U,n,m> &b) {
  Vector<typename Result<T,U>::Type,n> c;
  for (int k=0;k<n;k++) {
    c[k] = a[0]*b(k,0);
    for (int j=1;j<m;j++) c[k] += a[j]*b(k,j);
  }
  return c;
}


  // IO
template <typename T,int m,int n> inline
std::ostream &operator<<(std::ostream &os,const Matrix<T,m,n> &a) {
  for (int j=0;j<n;j++) {
    os << (j ? "]\n " : "\n[");
    for (int i=0;i<m;i++) os << (i?',':'[') << a(i,j);
  }
  os << "]]";
  return os;
}

template <typename T,int m,int n> inline
std::istream &operator>>(std::istream &is,Matrix<T,m,n> &a) {
  char ch;
  is >> ch;
  for (int j=0;j<n;j++) {
    is >> ch;
    for (int i=0;i<m;i++) is >> a(i,j) >> ch;
  }
  is >> ch;
  return is;
}


#ifdef Complex
#undef Complex
#endif

template <typename T>
class Complex {
public:
  Complex(void) {}
  Complex(const T &re) : re(re),im(0.0) {}
  Complex(const T &re,const T &im) : re(re),im(im) {}
  Complex(const T *const x) : re(x[0]),im(x[1]) {}
  static Complex FromRPhi(const T &r,const T &phi) {
    return Complex(r*cos(phi),r*sin(phi));
  }
  template <typename U> Complex(const Complex<U> &a) : re(a.re),im(a.im) {}
  template <typename U> const Complex &operator=(const Complex<U> &a)
    {re=a.re;im=a.im;return *this;}
  
    // access as 2-dimensional array: (TODO: find better hack)
  const T &operator[](int i) const {return (&re)[i];}
  T &operator[](int i) {return (&re)[i];}
  operator const T*(void) const {return &re;}
  operator T*(void) {return &re;}

    // special values:
  void setZero(void) {re=0;im=0;}
  bool nonZero(void) const {return ((re!=0) || (im!=0));}
  void setUnity(void) {re=1;im=0;}
  bool isUnity(void) const {return ((re==1) && (im==0));}

    // arithmetic
  template <typename U> const Complex &operator+=(const Complex<U> &a) {re+=a.re; im+=a.im; return *this;}
  template <typename U> const Complex &operator-=(const Complex<U> &a) {re-=a.re; im-=a.im; return *this;}
  template <typename U> const Complex &operator*=(const Complex<U> &b) {
    const T h = re*b.re - im*b.im;
    im        = im*b.re + re*b.im;
    re = h;
    return *this;
  }
  template <typename U> const Complex &operator/=(const Complex<U> &b) {return operator*=(b.inverse());}

  template <typename U> const Complex &operator*=(const U &x) {re*=x; im*=x; return *this;}
  template <typename U> const Complex &operator/=(const U &x) {re/=x; im/=x; return *this;}
  Complex operator-(void) const {return Complex(-re,-im);}
  const Complex &negate(void) {re=-re;im=-im;return *this;}
  Complex operator~(void) const {return Complex( re,-im);} // conjugation
  const Complex &conjugate(void) {im=-im;return *this;}
  Complex cross(void) const {return Complex(-im,re);} // multiplied by sqrt(-1)
  T length2(void) const {return re*re + im*im;}
  Complex inverse(void) const {return operator~() /= length2();}
  const Complex &invert(void) {
    const T h1(length2());
    re /= h1;
    im = -im/h1;
    return *this;
  }
public:
  T re,im;
};


template <typename T,typename U> inline
Complex<typename Result<T,U>::Type> operator+(const Complex<T> &v,
                                              const Complex<U> &w) {
  return Complex<typename Result<T,U>::Type>(v) += w;
}

template <typename T,typename U> inline
Complex<typename Result<T,U>::Type> operator-(const Complex<T> &v,
                                              const Complex<U> &w) {
  return Complex<typename Result<T,U>::Type>(v) -= w;
}

template <typename T,typename U> inline
Complex<typename MulResult<T,U>::Type> operator*(const Complex<T> &a,
                                                 const Complex<U> &b) {
  return Complex<typename Result<T,U>::Type>(a.re*b.re-a.im*b.im,a.im*b.re+a.re*b.im);
}

template <typename T,typename U> inline
Complex<typename Result<T,U>::Type> operator/(const Complex<T> &v,
                                              const Complex<U> &w) {
  return Complex<typename Result<T,U>::Type>(v) /= w;
}

template <typename T,typename U> inline
Complex<typename MulResult<T,U>::Type> operator*(const Complex<T> &a,
                                                 const U &x) {
  return Complex<typename MulResult<T,U>::Type>(a) *= x;
}

template <typename T,typename U> inline
Complex<typename MulResult<T,U>::Type> operator*(const U &x,
                                                 const Complex<T> &a) {
  return a*x;
}

template <typename T,typename U> inline
Complex<typename MulResult<T,U>::Type> operator/(const Complex<T> &a,
                                                 const U &x) {
  return Complex<typename MulResult<T,U>::Type>(a) /= x;
}

template <typename T> inline
Complex<T> sqr(const Complex<T> &a) {
  const T h = a.re*a.im;
  return Complex<T>(a.re*a.re-a.im*a.im,h+h);
}

template <typename T> inline
Complex<T> sqrt(const Complex<T> &a) {
    // rval.im >= 0
  if (abs(a.re) >= abs(a.im)) {
    const T re = sqrt(((T)(0.5))*(a.re + sqrt(a.length2())));
    return Complex<T>(re,((T)(0.5))*re/a.re);
  } else {
    const T im = copysign(sqrt(((T)(0.5))*(sqrt(a.length2()) - a.re)),a.im);
    return Complex<T>(((T)(0.5))*im/a.im,im);
  }
}

  // IO
template <typename T> inline
std::ostream &operator<<(std::ostream &os,const Complex<T> &a) {
  return os << '[' << a.re << ',' << a.im << ']';
}

template <typename T> inline
std::istream &operator>>(std::istream &is,Complex<T> &a) {
  char ch;
  if (!(is >> ch)) return is;
  if (ch != '[') {
    is.setstate(std::ios::failbit|std::ios::badbit);
    return is;
  }
  if (!(is >> a.re)) return is;
  if (!(is >> ch)) return is;
  if (ch != ',') {
    is.setstate(std::ios::failbit|std::ios::badbit);
    return is;
  }
  if (!(is >> a.im)) return is;
  if (!(is >> ch)) return is;
  if (ch != ']') is.setstate(std::ios::failbit|std::ios::badbit);
  return is;
}


template <typename T>
class Quaternion {
public:
    // type introspection
  typedef T Scal;
  typedef Vector<T,3> Vect;

    // construction
  Quaternion(void) {}
  Quaternion(const Scal &s) {scalar() = s;vector().setZero();}
  Quaternion(const Vect &v) {scalar() = 0;vector() = v;}
  Quaternion(const Scal &s,const Vect &v) {scalar() = s;vector() = v;}
  Quaternion(const Scal &s,const Scal *const v) {scalar() = s;vector() = v;}

  Quaternion(const Scal &x0,const Scal &x1,const Scal &x2,const Scal &x3)
    {scalar() = x0;vector(0) = x1;vector(1) = x2;vector(2) = x3;}
  Quaternion(const Scal *const x) {scalar() = x[0];vector() = x+1;}

  template <typename U> Quaternion(const Quaternion<U> &a)
    {scalar()=a.scalar();vector()=a.vector();}
  template <typename U> const Quaternion &operator=(const Quaternion<U> &a)
    {scalar()=a.scalar();vector()=a.vector();return *this;}

    // member access
  const Scal &scalar(void) const {return scal;}
  Scal &scalar(void) {return scal;}
  const Vect &vector(void) const {return vect;}
  Vect &vector(void) {return vect;}
  const Scal &vector(int i) const {return vect[i];}
  Scal &vector(int i) {return vect[i];}

    // access as 4-dimensional array: (TODO: find better hack)
  const T &operator[](int i) const {return (&scal)[i];}
  T &operator[](int i) {return (&scal)[i];}
  operator const T*(void) const {return &scal;}
  operator T*(void) {return &scal;}

    // special values:
  void setZero(void) {scalar()=0; vector().setZero();}
  bool nonZero(void) const {return ((scalar()!=0) || vector().nonZero());}
  void setUnity(void) {scalar()=1; vector().setZero();}
  bool isUnity(void) const {return ((scalar()==1) && (!vector().nonZero()));}

    // arithmetic
  template <typename U> const Quaternion &operator+=(const Quaternion<U> &a) {scalar()+=a.scalar(); vector()+=a.vector(); return *this;}
  template <typename U> const Quaternion &operator-=(const Quaternion<U> &a) {scalar()-=a.scalar(); vector()-=a.vector(); return *this;}
  template <typename U> const Quaternion &operator*=(const Quaternion<U> &a);

  template <typename U> const Quaternion &operator*=(const U &x) {scalar()*=x; vector()*=x; return *this;}
  template <typename U> const Quaternion &operator/=(const U &x) {scalar()/=x; vector()/=x; return *this;}
  Quaternion operator-(void) const {return Quaternion(-scalar(),-vector());}
  const Quaternion &negate(void) {scalar( )=-scalar( );vector().negate();return *this;}
  Quaternion operator~(void) const {return Quaternion( scalar(),-vector());} // conjugation
  const Quaternion &conjugate(void) {vector().negate();return *this;}
  Scal length2(void) const {return scalar( )*scalar( ) + vector().length2();}
  Quaternion inverse(void) const {return operator~() /= length2();}

  Matrix<T,3,3> rotationMatrix(void) const {
      // constructs the rotation matrix of the normalized Quaternion.
      // *this needs not to be already normalized.
    const Scal ww(scalar( )*scalar( ));
    const Scal wx(scalar( )*vector(0));
    const Scal wy(scalar( )*vector(1));
    const Scal wz(scalar( )*vector(2));
    const Scal xx(vector(0)*vector(0));
    const Scal xy(vector(0)*vector(1));
    const Scal xz(vector(0)*vector(2));
    const Scal yy(vector(1)*vector(1));
    const Scal yz(vector(1)*vector(2));
    const Scal zz(vector(2)*vector(2));
    const Scal f = 1/(ww+xx+yy+zz);
    const Scal f2 = f+f;
    return Matrix<T,3,3>(f*(ww+xx-yy-zz),f2*(-wz+xy),f2*(wy+xz),
                         f2*(wz+xy),f*(ww-xx+yy-zz),f2*(-wx+yz),
                         f2*(-wy+xz),f2*(wx+yz),f*(ww-xx-yy+zz));
  }
  template <typename U> Vector<typename Result<T,U>::Type,3>
                          rotateFromRight(const Vector<U,3> &u) const {
    const Scal ww(scalar( )*scalar( ));
    const Scal wx(scalar( )*vector(0));
    const Scal wy(scalar( )*vector(1));
    const Scal wz(scalar( )*vector(2));
    const Scal xx(vector(0)*vector(0));
    const Scal xy(vector(0)*vector(1));
    const Scal xz(vector(0)*vector(2));
    const Scal yy(vector(1)*vector(1));
    const Scal yz(vector(1)*vector(2));
    const Scal zz(vector(2)*vector(2));
    const Scal f = 1/(ww+xx+yy+zz);
    const Scal f2 = f+f;
    return Vector<typename Result<T,U>::Type,3>(
      (f*(ww+xx-yy-zz))*u[0] + (f2*(-wz+xy))*u[1] + (f2*(wy+xz))*u[2],
      (f2*(wz+xy))*u[0] + (f*(ww-xx+yy-zz))*u[1] + (f2*(-wx+yz))*u[2],
      (f2*(-wy+xz))*u[0] + (f2*(wx+yz))*u[1] + (f*(ww-xx-yy+zz))*u[2]);
  }
  template <typename U> Vector<typename Result<T,U>::Type,3>
                          rotateFromLeft(const Vector<U,3> &u) const {
    const Scal ww(scalar( )*scalar( ));
    const Scal wx(scalar( )*vector(0));
    const Scal wy(scalar( )*vector(1));
    const Scal wz(scalar( )*vector(2));
    const Scal xx(vector(0)*vector(0));
    const Scal xy(vector(0)*vector(1));
    const Scal xz(vector(0)*vector(2));
    const Scal yy(vector(1)*vector(1));
    const Scal yz(vector(1)*vector(2));
    const Scal zz(vector(2)*vector(2));
    const Scal f = 1/(ww+xx+yy+zz);
    const Scal f2 = f+f;
    return Vector<typename Result<T,U>::Type,3>(
      (f*(ww+xx-yy-zz))*u[0] + (f2*(wz+xy))*u[1] + (f2*(-wy+xz))*u[2],
      (f2*(-wz+xy))*u[0] + (f*(ww-xx+yy-zz))*u[1] + (f2*(wx+yz))*u[2],
      (f2*(wy+xz))*u[0] + (f2*(-wx+yz))*u[1] + (f*(ww-xx-yy+zz))*u[2]);
  }
  Quaternion(const Matrix<T,3,3> &rotation) {
      // hope no mistakes, did not try yet:
    const Scal h0 = rotation[0] + rotation[4] + rotation[8];
    const Scal h1 = rotation[0] - rotation[4] - rotation[8];
    const Scal h2 = rotation[4] - rotation[8] - rotation[0];
    const Scal h3 = rotation[8] - rotation[0] - rotation[4];
    if (h0 > h1) {
      if (h2 > h3) {
        if (h0 > h2) { // h0 is max
          goto h0_max;
        } else {       // h2 is max
          goto h2_max;
        }
      } else {
        if (h0 > h3) { // h0 is max
          h0_max:
          Scal h = sqrt(1 + h0);
          scalar( ) = 0.5 * h;
          h         = 0.5 / h;
          vector(0) = h * (rotation[7] - rotation[5]);
          vector(1) = h * (rotation[2] - rotation[6]);
          vector(2) = h * (rotation[3] - rotation[1]);
        } else {       // h3 is max
          goto h3_max;
        }
      }
    } else {
      if (h2 > h3) {
        if (h1 > h2) { // h1 is max
          goto h1_max;
        } else {       // h2 is max
          h2_max:
          Scal h  = sqrt(1 + h2);
          vector(1) = 0.5 * h;
          h         = 0.5 / h;
          vector(2) = h * (rotation[7] + rotation[5]);
          scalar( ) = h * (rotation[2] - rotation[6]);
          vector(0) = h * (rotation[3] + rotation[1]);
        }
      } else {
        if (h1 > h3) { // h1 is max
          h1_max:
          Scal h  = sqrt(1 + h1);
          vector(0) = 0.5 * h;
          h         = 0.5 / h;
          scalar( ) = h * (rotation[7] - rotation[5]);
          vector(2) = h * (rotation[2] + rotation[6]);
          vector(1) = h * (rotation[3] + rotation[1]);
        } else {       // h3 is max
          h3_max:
          Scal h  = sqrt(1 + h3);
          vector(2) = 0.5 * h;
          h         = 0.5 / h;
          vector(1) = h * (rotation[7] + rotation[5]);
          vector(0) = h * (rotation[2] + rotation[6]);
          scalar( ) = h * (rotation[3] - rotation[1]);
        }
      }
    }
  }
private:
  Scal scal;
  Vect vect;
};


template <typename T,typename U> inline
Quaternion<typename Result<T,U>::Type> operator+(const Quaternion<T> &v,
                                                 const Quaternion<U> &w) {
  return Quaternion<typename Result<T,U>::Type>(v) += w;
}

template <typename T,typename U> inline
Quaternion<typename Result<T,U>::Type> operator-(const Quaternion<T> &v,
                                                 const Quaternion<U> &w) {
  return Quaternion<typename Result<T,U>::Type>(v) -= w;
}

template <typename T,typename U> inline
Quaternion<typename MulResult<T,U>::Type> operator*(const Quaternion<T> &a,
                                                    const Quaternion<U> &b) {
  return Quaternion<typename Result<T,U>::Type>(
           a.scalar()*b.scalar( ) - a.vector(0)*b.vector(0) - a.vector(1)*b.vector(1) - a.vector(2)*b.vector(2),
           a.scalar()*b.vector(0) + a.vector(0)*b.scalar( ) + a.vector(1)*b.vector(2) - a.vector(2)*b.vector(1),
           a.scalar()*b.vector(1) - a.vector(0)*b.vector(2) + a.vector(1)*b.scalar( ) + a.vector(2)*b.vector(0),
           a.scalar()*b.vector(2) + a.vector(0)*b.vector(1) - a.vector(1)*b.vector(0) + a.vector(2)*b.scalar( ));
}

template<class T> template<class U> 
const Quaternion<T> &Quaternion<T>::operator*=(const Quaternion<U> &a) {
  Quaternion<T> r = (*this)*a;
  operator=(r);
  return *this;
}

template <typename T,typename U> inline
Quaternion<typename MulResult<T,U>::Type> operator*(const Quaternion<T> &a,
                                                    const U &x) {
  return Quaternion<typename MulResult<T,U>::Type>(a) *= x;
}

template <typename T,typename U> inline
Quaternion<typename MulResult<T,U>::Type> operator*(const U &x,
                                                    const Quaternion<T> &a) {
  return a*x;
}

template <typename T,typename U> inline
Quaternion<typename MulResult<T,U>::Type> operator/(const Quaternion<T> &a,
                                                    const U &x) {
  return Quaternion<typename MulResult<T,U>::Type>(a) /= x;
}

template <typename T,typename U> inline
Quaternion<typename MulResult<T,U>::Type> operator*(const Quaternion<T> &a,
                                                    const Vector<U,3> &b) {
  return Quaternion<typename Result<T,U>::Type>(
         - a.vector(0)*b[0] - a.vector(1)*b[1] - a.vector(2)*b[2],
           a.scalar( )*b[0] - a.vector(2)*b[1] + a.vector(1)*b[2],
           a.vector(2)*b[0] + a.scalar( )*b[1] - a.vector(0)*b[2],
         - a.vector(1)*b[0] + a.vector(0)*b[1] + a.scalar( )*b[2]);
}

template <typename T,typename U> inline
Quaternion<typename MulResult<T,U>::Type> operator*(const Vector<T,3> &a,
                                                    const Quaternion<U> &b) {
  return Quaternion<typename Result<T,U>::Type>(
           - a[0]*b.vector(0) - a[1]*b.vector(1) - a[2]*b.vector(2),
             a[0]*b.scalar( ) + a[1]*b.vector(2) - a[2]*b.vector(1),
           - a[0]*b.vector(2) + a[1]*b.scalar( ) + a[2]*b.vector(0),
             a[0]*b.vector(1) - a[1]*b.vector(0) + a[2]*b.scalar( ));
}





  // IO
template <typename T> inline
std::ostream &operator<<(std::ostream &os,const Quaternion<T> &a) {
  return os << '[' << a.scalar() << ',' << a.vector() << ']';
}

template <typename T> inline
std::istream &operator>>(std::istream &is,Quaternion<T> &a) {
  char ch;
  is >> ch >> a.scalar() >> ch >> a.vector() >> ch;
  return is;
}


#endif
