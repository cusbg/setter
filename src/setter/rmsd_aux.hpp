#ifndef RMSD_AUX_HPP_
#define RMSD_AUX_HPP_

#include "rmsd.h"

#include <iostream>
#include <fstream>

/**
 * The cRMSD3DCoord class represents a 3-dimensional vector.
 */
class cRMSD3DCoord
{
private:
    double data[3];
    
public:
    cRMSD3DCoord()
    {
        for(int i = 0; i < 3; i++)
            data[i] = 0.0;
    }
    
    cRMSD3DCoord(double value[3])
    {
        for(int i = 0; i < 3; i++)
            data[i] = value[i];
    }
    
    cRMSD3DCoord(double x, double y, double z)
    {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    
    cRMSD3DCoord(const cRMSD3DCoord &v)
    {
        for(int i = 0; i < 3; i++)
            data[i] = v.data[i];
    }
    
    double &
    operator[](int index)
    {
        return data[index];
    }

};



inline cRMSD3DCoord
operator+(cRMSD3DCoord a, cRMSD3DCoord b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        v[i] = a[i] + b[i];
    
    return v;
}



inline cRMSD3DCoord
operator-(cRMSD3DCoord a, cRMSD3DCoord b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        v[i] = a[i] - b[i];
    
    return v;
}



inline double
operator*(cRMSD3DCoord a, cRMSD3DCoord b)
{
    double d = 0.0;
    
    for(int i = 0; i < 3; i++)
        d += a[i] * b[i];
    
    return d;
}



inline cRMSD3DCoord
operator*(double a, cRMSD3DCoord b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        v[i] = a * b[i];
    
    return v;
}



inline cRMSD3DCoord
operator*(cRMSD3DCoord a, double b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        v[i] = a[i] * b;
    
    return v;
}



inline cRMSD3DCoord
operator/(cRMSD3DCoord a, double b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        v[i] = a[i] / b;
    
    return v;
}

inline std::ostream &
operator<<(std::ostream &stream, cRMSD3DCoord a)
{
    stream << "{ ";
    
    for(int i = 0; i < 3; i++)
        stream << a[i] << " ";
    
    stream << "}";
    return stream;
}



/**
 * Computes the square of the norm of the vector.
 */
inline double
norm2(cRMSD3DCoord a)
{
    return a*a;
}

// *********** MATRIX *********************

/**
 * The cRMSDMatrix class represents a 3-by-3 matrix. 
 */
class cRMSDMatrix
{
private:
    cRMSD3DCoord data[3];

public:
    cRMSDMatrix()
    {
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                data[i][j] = 0.0;
    }
    
    cRMSDMatrix(double value[3][3])
    {
        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                data[i][j] = value[i][j];
    }
    
    cRMSDMatrix(const cRMSDMatrix &v)
    {
        for(int i = 0; i < 3; i++)
            data[i] = v.data[i];
    }
    
    cRMSD3DCoord &
    operator[](int index)
    {
        return data[index];
    }
};



inline cRMSDMatrix
operator+(cRMSDMatrix a, cRMSDMatrix b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        m[i] = a[i] + b[i];
    
    return m;
}



inline cRMSDMatrix
operator-(cRMSDMatrix a, cRMSDMatrix b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        m[i] = a[i] - b[i];
    
    return m;
}

inline double
operator*(cRMSDMatrix a, cRMSDMatrix b)
{
    double product = 0;
    
    for(int i = 0; i < 3; i++)
        product += a[i]*b[i];
    
    return product;
}



inline cRMSDMatrix
operator*(double a, cRMSDMatrix b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        m[i] = a * b[i];
    
    return m;
}



inline cRMSDMatrix
operator*(cRMSDMatrix a, double b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        m[i] = a[i] * b;
    
    return m;
}



inline cRMSD3DCoord
operator*(cRMSDMatrix a, cRMSD3DCoord b)
{
    cRMSD3DCoord v;
    
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            v[i] += a[i][j] * b[j];
    
    return v;
}



/**
 * Computes the vector direct product of vectors.
 */
inline cRMSDMatrix
operator^(cRMSD3DCoord a, cRMSD3DCoord b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            m[i][j] = a[i] * b[j];
    
    return m;
}



inline cRMSDMatrix
operator/(cRMSDMatrix a, double b)
{
    cRMSDMatrix m;
    
    for(int i = 0; i < 3; i++)
        m[i] = a[i] / b;
    
    return m;
}


inline double
norm2(cRMSDMatrix a)
{
	return norm2(a[0] + a[1] + a[2])/3;
}

// **************** END MATRIX *****************


/**
 * The cRMSDStruct class represents a protein structure.
 */
class cRMSDStruct
{
public:
	std::string name;           // name of the protein (need for generate_vmd)
    
    unsigned int length;                 // number of resids (max length)
	unsigned int lengthTemp;             // number of resids (lengthTemp <= length)
    cRMSD3DCoord *coord;              // coordinates of resids
    std::string aa;
	std::string sse;                  // type of secondary structure (A,B,T)
    
    int resid_offset;           // number of first resid in PDB (need for generate_vmd)
    
    int family;
    int superfamily;
    int fold;
    int klass;
    
    cRMSDStruct()
    {
		//name = new char[100]; 
		length = 0;
		coord = NULL;
    }
    
    /**
     * Creates a new transformed protein structure.
     * 
     * \param p     template protein
     * \param u     rotation
     * \param t     transposition
     */
    cRMSDStruct(cRMSDStruct *p, const cRMSDMatrix &u, const cRMSD3DCoord &t)
    {
		//name = new char[100]; 

		SetName(p->name);
		//strcpy(name, p->name); 
        length = p->length;
		lengthTemp = p->lengthTemp;
        aa = p->aa;
		sse = p->sse;
        resid_offset = p->superfamily;
        superfamily = p->superfamily;
        fold = p->fold;
        klass = p->klass;
        
        coord = new cRMSD3DCoord[length];
        
        for(unsigned int i = 0; i < length; i++)
            coord[i] = u*p->coord[i] + t;
    }	

	void Reinit()
	{
		if (coord != NULL)
		{
			delete[] coord;
			coord = NULL;
			length = 0;
			name = "";
		}
	}
    
    ~cRMSDStruct()
    {
		if (coord != NULL)
		{
			delete[] coord;			
		}
    }

	void SetName(std::string _name)
	{
		name = _name;
	}
    
    /*
    int
    operator()(cRMSD3DCoord &v)
    {
        int res = -1;
        double best = 0;
        
        for(int i = 0; i < length; i++)
        {
            double d = norm2(coord[i] - v);
            
            if(d < best || res == -1)
            {
                res = i;
                best = d;
            }
        }
        
        return res;
    }
    */
};


// ****************** ALIGN *********************

#define NOWEIGHT

/**
 * Weight of the distance of the resid pair.
 * 
 * \param qprot         query protein
 * \param dprot         database protein
 * \param qres          resid in the query protein
 * \param dres          resid in the database protein
 */
inline double
pair_weight(cRMSDStruct *qprot, cRMSDStruct *dprot, int qres, int dres) {
    return 1;
}



/**
 * Alignment item
 */
typedef struct
{
    int x;      // resid number in query protein (the first has 0)
    int y;      // resid number in database protein (the first has 0)
#ifndef NOWEIGHT
    double w;   // weight
#endif /*NOWEIGHT*/    
}
align_t;



inline bool
operator!=(align_t a, align_t b)
{
    if(a.x != b.x || a.y != b.y)
        return true;
    else
        return false;
}



/**
 * Alignment
 */
class cRMSDAlign
{
    align_t *al;    // data
    int len;        // length
    int *counter;   // data reference counter
    
public:
    
    /**
     * Create an empty alignment with the given length. 
     */
    cRMSDAlign(int length=0) : al(0), counter(0)
    {
        len = length;
        
        if(len)
        {
            al = new align_t[len];
            counter = new int(1);
        }
    }
    
    
    /**
     * Create a shallow copy.
     */
    cRMSDAlign(const cRMSDAlign &align)
    {
        al      = align.al;
        len     = align.len;
        counter = align.counter;
        
        if(counter)
            (*counter)++;

    }
    
    
    /**
     * Create an alignment from a DP path. 
     */
    cRMSDAlign(int plen, int (*path)[2]);
    
    
    /**
     * Destructor.
     */
    ~cRMSDAlign()
    {
        if(counter)
        {
            (*counter)--;
            
            if(! *counter)
            {
                delete counter;
                delete[] al;
            }
        }
    }
    
    
    /**
     * Assignment operator (shallow copy)
     */
    cRMSDAlign &operator=(const cRMSDAlign &align)
    {
        if(counter && al == align.al)
            return *this;
        
        if(counter)
        {
            (*counter)--;
            
            if(! *counter)
            {
                delete counter;
                delete[] al;
            }
        }
        
        al      = align.al;
        len     = align.len;
        counter = align.counter;
        
        if(counter)
            (*counter)++;
        
        return *this;
    }
    
    
    /**
     * Fill the alignment according to a part of other alignment.
     */
    void
    set(cRMSDAlign &base, int from, int to)
    {
        len = to - from;
        
        for(int i = from; i < to; i++)
        {
            al[i-from].x = base.al[i].x;
            al[i-from].y = base.al[i].y;
#ifndef NOWEIGHT
            al[i-from].w = base.al[i].w;
#endif /*NOWEIGHT*/    

        }
    }
    
    
    /**
     * Add to the alignment a part of other alignment.
     */
    void
    add(cRMSDAlign &base, int from, int to)
    {
        int offset = len;
        
        len += to - from;
        
        for(int i = from; i < to; i++)
        {
            al[offset + i-from].x = base.al[i].x;
            al[offset + i-from].y = base.al[i].y;
#ifndef NOWEIGHT
            al[offset + i-from].w = base.al[i].w;
#endif /*NOWEIGHT*/    
        }
    }
    
    
    /**
     * Returns the length of the alignment.
     */
    int
    length()
    {
        return len;
    }
    
    
    /**
     * Set length of the alignment.
     */
    void
    setlength(int len)
    {
        this->len = len;
    }
    
    
    align_t &
    operator[](int i)
    {
        return al[i];
    }
    
    
    bool
    operator==(cRMSDAlign &align)
    {
        return al == align.al;
    }
    
    
    cRMSDAlign
    copy()
    {
        cRMSDAlign align(len);
        align.set(*this, 0, len);
        return align;
    }
    
    
    /**
     * Return an alignment after DP phase.
     */
    cRMSDAlign dp(cRMSDStruct *cprot, cRMSDStruct *dprot, int dpwidth, double d0);
};   

// ****************** END ALIGN *********************

// ****************** RMSD *********************

double rmsd(cRMSDStruct &s1, cRMSDStruct &s2, cRMSDAlign &align, cRMSDMatrix &rot, cRMSD3DCoord &trans)
{
	double structRef[10000][3];
	double structMov[10000][3];
	/*structRef = new (double[3])[align.length()];
	structMov = new double*[align.length()];
	for (int i = 0; i < align.length(); i++)
	{
		structRef[i] = new double[3];
		structMov[i] = new double[3];
	}*/

	for (int i = 0; i < align.length(); i++)
	{ 
		for (int j = 0; j < 3; j++)
		{
			structRef[i][j] = s1.coord[align[i].x][j]; 
			structMov[i][j] = s2.coord[align[i].y][j];
		}
	}

	double mov_com[3];
	double mov_to_ref[3];
	double U[3][3];
	double rmsd_dist;

	calculate_rotation_rmsd(structRef, structMov, align.length(), mov_com, mov_to_ref, U, &rmsd_dist);	

	rot = cRMSDMatrix(U);
	trans = cRMSD3DCoord(mov_to_ref);

	/*for (int i = 0; i < align.length(); i++)
	{
		delete [] structRef[i];
		delete [] structMov[i];
	}
	delete [] structRef;
	delete [] structMov;*/

	return rmsd_dist;
}

#endif /* RMSD_AUX_HPP_ */