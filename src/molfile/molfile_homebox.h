#include <math.h>

static inline
double molfile_dotprod(const double* x, const double* y) {
    return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

static inline
void molfile_unitcell_to_pdb( const double* m_box,
                              double* a, double* b, double* c,
                              double* alpha, double* beta, double* gamma) {

    double A[3] = { m_box[0], m_box[1], m_box[2] };
    double B[3] = { m_box[3], m_box[4], m_box[5] };
    double C[3] = { m_box[6], m_box[7], m_box[8] };

    *a = *b = *c = 1;
    *alpha = *beta = *gamma = 90;

    /* store lengths */
    *a = sqrt(molfile_dotprod(A,A));
    *b = sqrt(molfile_dotprod(B,B));
    *c = sqrt(molfile_dotprod(C,C));

    if (*a && *b && *c) {
        /* compute angles */
        double cosAB = molfile_dotprod(A,B)/(a[0] * b[0]);
        double cosAC = molfile_dotprod(A,C)/(a[0] * c[0]);
        double cosBC = molfile_dotprod(B,C)/(b[0] * c[0]);

        // clamp
        if (cosAB > 1.0) cosAB = 1.0; else if (cosAB < -1.0) cosAB = -1.0;
        if (cosAC > 1.0) cosAC = 1.0; else if (cosAC < -1.0) cosAC = -1.0;
        if (cosBC > 1.0) cosBC = 1.0; else if (cosBC < -1.0) cosBC = -1.0;

        /* convert to angles using asin to avoid nasty rounding when we are */
        /* close to 90 degree angles. */
        *alpha = 90.0 - asin(cosBC) * 90.0 / M_PI_2; /* cosBC */
        *beta  = 90.0 - asin(cosAC) * 90.0 / M_PI_2; /* cosAC */
        *gamma = 90.0 - asin(cosAB) * 90.0 / M_PI_2; /* cosAB */
    }
}

static inline
void molfile_unitcell_from_pdb( double* m_box,
                                double a, double b, double c,
                                double alpha, double beta, double gamma) {

    double A[3], B[3], C[3];

    // Convert VMD's unit cell information
    double cosBC = sin( ((90 - alpha ) / 180) * M_PI );
    double cosAC = sin( ((90 - beta  ) / 180) * M_PI );
    double cosAB = sin( ((90 - gamma ) / 180) * M_PI );
    double sinAB = cos( ((90 - gamma ) / 180) * M_PI );

    double Ax = a;
    double Ay = 0;
    double Az = 0;
    double Bx = b * cosAB;
    double By = b * sinAB;
    double Bz = 0;
    double Cx,Cy,Cz;
    if (sinAB != 0) {
        Cx = cosAC;
        Cy = (cosBC - cosAC*cosAB) / sinAB;
        Cz = sqrt(1-Cx*Cx-Cy*Cy);
        Cx *= c;
        Cy *= c;
        Cz *= c;
    } else {
        Cx=Cy=Cz=0;
    }
    A[0] = Ax; A[1] = Ay; A[2] = Az;
    B[0] = Bx; B[1] = By; B[2] = Bz;
    C[0] = Cx; C[1] = Cy; C[2] = Cz;

    /* put vectors in rows of box */
    m_box[0] = A[0]; m_box[1] = A[1]; m_box[2] = A[2];
    m_box[3] = B[0]; m_box[4] = B[1]; m_box[5] = B[2];
    m_box[6] = C[0]; m_box[7] = C[1]; m_box[8] = C[2];
}
