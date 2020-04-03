
#ifndef IMAGE3D_H
#define IMAGE3D_H



//RFBA: Image3D Class for 3D Image Processing, Added@20180306


#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <memory>
#include <algorithm>

//Region3D
class Ref3D {
public:
	int x;
	int y;
	int z;

	inline Ref3D() {
		x = 0;
		y = 0;
		z = 0;
	};

	Ref3D(int xin, int yin, int zin) {
		x = xin;
		y = yin;
		z = zin;
	}

	Ref3D(const Ref3D& size) {
		x = size.x;
		y = size.y;
		z = size.z;
	}

	Ref3D(Ref3D& size) {
		x = size.x;
		y = size.y;
		z = size.z;
	}



	inline bool Ref3D::operator<(const Ref3D & other) const
	{
		return z < other.z || (z == other.z && y < other.y) || (z == other.z && y == other.y && x < other.x);
	}

	inline bool Ref3D::operator>(const Ref3D & other) const
	{
		return z > other.z || (z == other.z && y > other.y) || (z == other.z && y == other.y && x > other.x);
	}

	inline Ref3D Ref3D::operator+(const Ref3D rhs) const
	{
		Ref3D v(x + rhs.x, y + rhs.y, z + rhs.z);
		return v;
	}

	inline bool Ref3D::operator ==(const Ref3D& ref) const
	{
		return (x == ref.x && y == ref.y && z == ref.z);
	}

	inline bool Ref3D::operator !=(const Ref3D& ref) const
	{
		return (x != ref.x || y != ref.y || z != ref.z);
	}


	//Fuctions
	inline unsigned int Ref3D::mag_squared() const
	{
		typedef unsigned int uint;
		return uint(x*x) + uint(y*y) + uint(z*z);
	}


	inline int Ref3D::area() {
		return(x*y*z);
	}

};



//Image 3D
template<class T> class Image3D
{
protected:
	T* my_data;       ///< The raw image data
	Ref3D my_size; ///< The size of the image


public:

	Image3D(const Ref3D& size) {
		my_size = Ref3D(size);
		my_data = (T*)malloc(my_size.area() * sizeof(T));
		for (int i = 0; i < my_size.area(); i++)
			my_data[i] = 0;
	}

	Image3D(const Ref3D& size, const T val) {
		my_size = Ref3D(size);
		my_data = (T*)malloc(my_size.area() * sizeof(T));
		for (int i = 0; i < my_size.area(); i++)
			my_data[i] = val;
	}

	Image3D(const Ref3D& size, const T* input) {
		my_size = Ref3D(size);
		my_data = (T*)malloc(my_size.area() * sizeof(T));
		for (int i = 0; i < my_size.area(); i++)
			my_data[i] = input[i];
	}


	inline Ref3D size() const
	{
		return my_size;
	}


	/// Access a pixel from the image. Bounds checking is only performed if the library is compiled
	/// with <code>-D CVD_IMAGE_DEBUG</code>, in which case an ImageError::AccessOutsideImage exception is 
	/// thrown.
	inline T& operator[](const Ref3D& pos)
	{
		return (my_data[pos.z*my_size.x*my_size.y + pos.y*my_size.x + pos.x]);
	}
	/// Access a pixel from the image. Bounds checking is only performed if the library is compiled
	/// with <code>-D CVD_IMAGE_DEBUG</code>, in which case an ImageError::AccessOutsideImage exception is 
	/// thrown.
	inline const T& operator[](const Ref3D& pos) const
	{
		return (my_data[pos.z*my_size.x*my_size.y + pos.y*my_size.x + pos.x]);
	}
	inline const T& operator[](const int pos) const
	{
		return (my_data[pos]);
	}
	/// Returns the raw image data
	inline const T* data() const
	{
		return(my_data);
	}

	/// Returns the raw image data
	inline T* data()
	{
		return(my_data);
	}

	inline  Image3D<T> copy_from_me() const {
		Image3D<T> img(my_size);
		for (int i = 0; i < my_size.x * my_size.y * my_size.z; i++)
			img.data()[i] = my_data[i];
		return(img);
	}


};



#endif