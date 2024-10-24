/**
 *
 * Original author: Duc Nguyen
 */

#include "modules/ddc/ddc.h"

namespace chflow {

int field2vector_size(const FlowField& u, const FlowField& temp, const FlowField& salt) {
    int Kx = u.kxmaxDealiased();
    int Kz = u.kzmaxDealiased();
    int Ny = u.Ny();
    // FIXME: Determine array size
    // The original formula was
    //     int N = 4* ( Kx+2*Kx*Kz+Kz ) * ( Ny-3 ) +2* ( Ny-2 );
    // but I've not been able to twist my head enough to adapt it do distributed FlowFields.
    // Since it doesn't take that much time, we now perform the loops twice, once empty to
    // determine cfarray sizes and once with the actual data copying. // Tobias
    int N = 0;
    if (u.taskid() == u.task_coeff(0, 0))
        N += 2 * (Ny - 2);
    for (int kx = 1; kx <= Kx; ++kx)
        if (u.taskid() == u.task_coeff(u.mx(kx), 0))
            N += 2 * (Ny - 2) + 2 * (Ny - 4);
    for (int kz = 1; kz <= Kz; ++kz)
        if (u.taskid() == u.task_coeff(0, u.mz(kz)))
            N += 2 * (Ny - 2) + 2 * (Ny - 4);
    for (int kx = -Kx; kx <= Kx; ++kx) {
        if (kx == 0)
            continue;
        int mx = u.mx(kx);
        for (int kz = 1; kz <= Kz; ++kz) {
            int mz = u.mz(kz);
            if (u.taskid() == u.task_coeff(mx, mz)) {
                N += 2 * (Ny - 2) + 2 * (Ny - 4);
            }
        }
    }

    if (temp.taskid() == temp.task_coeff(0, 0))
        N += Ny;
    for (int kx = 1; kx <= Kx; ++kx)
        if (temp.taskid() == temp.task_coeff(temp.mx(kx), 0))
            N += 2 * Ny;
    for (int kz = 1; kz <= Kz; ++kz)
        if (temp.taskid() == temp.task_coeff(0, temp.mz(kz)))
            N += 2 * Ny;
    for (int kx = -Kx; kx <= Kx; kx++) {
        if (kx == 0)
            continue;
        for (int kz = 1; kz <= Kz; kz++) {
            if (temp.taskid() == temp.task_coeff(temp.mx(kx), temp.mz(kz))) {
                N += 2 * Ny;
            }
        }
    }

    if (salt.taskid() == salt.task_coeff(0, 0))
        N += Ny;
    for (int kx = 1; kx <= Kx; ++kx)
        if (salt.taskid() == salt.task_coeff(salt.mx(kx), 0))
            N += 2 * Ny;
    for (int kz = 1; kz <= Kz; ++kz)
        if (salt.taskid() == salt.task_coeff(0, salt.mz(kz)))
            N += 2 * Ny;
    for (int kx = -Kx; kx <= Kx; kx++) {
        if (kx == 0)
            continue;
        for (int kz = 1; kz <= Kz; kz++) {
            if (salt.taskid() == salt.task_coeff(salt.mx(kx), salt.mz(kz))) {
                N += 2 * Ny;
            }
        }
    }

    return N;
}

/** \brief Turn the three flowfields for velocity, temperature, and salinity into one Eigen vector
 * \param[in] u velocity field
 * \param[in] temp temperature field
 * \param[in] salt salinity field
 * \param[in] a vector for the linear algebra
 *
 * The vectorization of u is analog to the field2vector, temparture and salinity are piped entirely
 * into the vector (a single independent dimensions)
 */
void field2vector(const FlowField& u, const FlowField& temp, const FlowField& salt, Eigen::VectorXd& a) {
    Eigen::VectorXd b;
    assert(temp.xzstate() == Spectral && temp.ystate() == Spectral);
    assert(salt.xzstate() == Spectral && salt.ystate() == Spectral);
    field2vector(u, b);
    int Kx = temp.kxmaxDealiased();
    int Kz = temp.kzmaxDealiased();
    int Ny = temp.Ny();

    int n = field2vector_size(u, temp, salt);  // b.size() +6*Kx*Kz*Ny;

    if (a.size() < n)
        a.resize(n, true);
    setToZero(a);
    int pos = b.size();
    a.topRows(pos) = b;

    // (0,:,0)
    for (int ny = 0; ny < Ny; ++ny)
        if (temp.taskid() == temp.task_coeff(0, 0))
            a(pos++) = Re(temp.cmplx(0, ny, 0, 0));

    for (int kx = 1; kx <= Kx; ++kx) {
        int mx = temp.mx(kx);
        if (temp.taskid() == temp.task_coeff(mx, 0)) {
            for (int ny = 0; ny < Ny; ++ny) {
                a(pos++) = Re(temp.cmplx(mx, ny, 0, 0));
                a(pos++) = Im(temp.cmplx(mx, ny, 0, 0));
            }
        }
    }
    for (int kz = 1; kz <= Kz; ++kz) {
        int mz = temp.mz(kz);
        if (temp.taskid() == temp.task_coeff(0, mz)) {
            for (int ny = 0; ny < Ny; ++ny) {
                a(pos++) = Re(temp.cmplx(0, ny, mz, 0));
                a(pos++) = Im(temp.cmplx(0, ny, mz, 0));
            }
        }
    }
    for (int kx = -Kx; kx <= Kx; kx++) {
        if (kx == 0)
            continue;
        int mx = temp.mx(kx);
        for (int kz = 1; kz <= Kz; kz++) {
            int mz = temp.mz(kz);
            if (u.taskid() == u.task_coeff(mx, mz)) {
                for (int ny = 0; ny < Ny; ny++) {
                    a(pos++) = Re(temp.cmplx(mx, ny, mz, 0));
                    a(pos++) = Im(temp.cmplx(mx, ny, mz, 0));
                }
            }
        }
    }
    // <<----- missing salt
}

/** \brief Turn  one Eigen vector into the three flowfields for velocity, temperature, and salinity
 * \param[in] a vector for the linear algebra
 * \param[in] u velocity field
 * \param[in] temp temperature field
 * \param[in] salt temperature field
 *
 * The vectorization of u is analog to the field2vector, temperature and salinity are piped entirely
 * into the vector (a single independent dimension)
 */
void vector2field(const Eigen::VectorXd& a, FlowField& u, FlowField& temp, FlowField& salt) {
    assert(temp.xzstate() == Spectral && temp.ystate() == Spectral);
    assert(salt.xzstate() == Spectral && salt.ystate() == Spectral);
    temp.setToZero();
    Eigen::VectorXd b;
    int N = field2vector_size(u);
    int Kx = temp.kxmaxDealiased();
    int Kz = temp.kzmaxDealiased();
    int Ny = temp.Ny();
    b = a.topRows(N);
    vector2field(b, u);

    int pos = N;
    double reval, imval;
    Complex val;

    if (temp.taskid() == temp.task_coeff(0, 0))
        for (int ny = 0; ny < Ny; ++ny)
            temp.cmplx(0, ny, 0, 0) = Complex(a(pos++), 0);

    for (int kx = 1; kx <= Kx; ++kx) {
        int mx = temp.mx(kx);
        if (temp.taskid() == temp.task_coeff(mx, 0)) {
            for (int ny = 0; ny < Ny; ++ny) {
                reval = a(pos++);
                imval = a(pos++);
                temp.cmplx(mx, ny, 0, 0) = Complex(reval, imval);
            }
        }

        // ------------------------------------------------------
        // Now copy conjugates of u(kx,ny,0,i) to u(-kx,ny,0,i). These are
        // redundant modes stored only for the convenience of FFTW.
        int mxm = temp.mx(-kx);
        int send_id = temp.task_coeff(mx, 0);
        int rec_id = temp.task_coeff(mxm, 0);
        for (int ny = 0; ny < Ny; ++ny) {
            if (temp.taskid() == send_id && send_id == rec_id) {  // all is on the same process -> just copy
                temp.cmplx(mxm, ny, 0, 0) = conj(temp.cmplx(mx, ny, 0, 0));
            }
#ifdef HAVE_MPI     // send_id != rec_id requires multiple processes
            else {  // Transfer the conjugates via MPI
                if (temp.taskid() == send_id) {
                    Complex tmp0 = conj(temp.cmplx(mx, ny, 0, 0));
                    MPI_Send(&tmp0, 1, MPI_DOUBLE_COMPLEX, rec_id, 0, MPI_COMM_WORLD);
                }
                if (u.taskid() == rec_id) {
                    Complex tmp0;
                    MPI_Status status;
                    MPI_Recv(&tmp0, 1, MPI_DOUBLE_COMPLEX, send_id, 0, MPI_COMM_WORLD, &status);
                    temp.cmplx(mxm, ny, 0, 0) = tmp0;
                }
            }
#endif
        }
    }
    for (int kz = 1; kz <= Kz; ++kz) {
        int mz = temp.mz(kz);
        if (temp.taskid() == temp.task_coeff(0, mz)) {
            for (int ny = 0; ny < Ny; ++ny) {
                reval = a(pos++);
                imval = a(pos++);
                temp.cmplx(0, ny, mz, 0) = Complex(reval, imval);
            }
        }
    }
    for (int kx = -Kx; kx <= Kx; kx++) {
        if (kx == 0)
            continue;
        int mx = temp.mx(kx);
        for (int kz = 1; kz <= Kz; kz++) {
            int mz = temp.mz(kz);
            if (u.taskid() == u.task_coeff(mx, mz)) {
                for (int ny = 0; ny < Ny; ny++) {
                    reval = a(pos++);
                    imval = a(pos++);
                    val = Complex(reval, imval);
                    temp.cmplx(mx, ny, mz, 0) = val;
                }
            }
        }
    }
    temp.setPadded(true);


    // <<----- missing salt
}

// DDC::DDC()
//   :
//   DNS(){
//
// }

DDC::DDC(const std::vector<FlowField>& fields, const DDCFlags& flags)
    :  // base class constructor with no arguments is called automatically (see DNS::DNS())
      main_dde_(0),
      init_dde_(0) {
    main_dde_ = newDDE(fields, flags);
    // creates DNSAlgo with ptr of "nse"-daughter type "dde"
    main_algorithm_ = newAlgorithm(fields, main_dde_, flags);
    if (!main_algorithm_->full() && flags.initstepping != flags.timestepping) {
        DDCFlags initflags = flags;
        initflags.timestepping = flags.initstepping;
        initflags.dt = flags.dt;
        init_dde_ = newDDE(fields, flags);

        // creates DNSAlgo with ptr of "nse"-daughter type "dde"
        init_algorithm_ = newAlgorithm(fields, init_dde_, initflags);
        // Safety check
        if (init_algorithm_->Ninitsteps() != 0)
            std::cerr << "DDC::DDC(fields, flags) :\n"
                      << flags.initstepping << " can't initialize " << flags.timestepping
                      << " since it needs initialization itself.\n";
    }
    
}

DDC::~DDC() {}



const ChebyCoeff& DDC::Ubase() const {
    if (main_dde_)
        return main_dde_->Ubase();
    else if (init_dde_)
        return init_dde_->Ubase();
    else {
        std::cerr << "Error in DDC::Ubase(): Ubase is currently undefined" << std::endl;
        exit(1);
        return init_dde_->Ubase();  // to make compiler happy
    }
}
const ChebyCoeff& DDC::Wbase() const {
    if (main_dde_)
        return main_dde_->Wbase();
    else if (init_dde_)
        return init_dde_->Wbase();
    else {
        std::cerr << "Error in DDC::Wbase(): Wbase is currently undefined" << std::endl;
        exit(1);
        return init_dde_->Wbase();  // to make compiler happy
    }
}
const ChebyCoeff& DDC::Tbase() const {
    if (main_dde_)
        return main_dde_->Tbase();
    else if (init_dde_)
        return init_dde_->Tbase();
    else {
        std::cerr << "Error in DDC::Tbase(): Tbase is currently undefined" << std::endl;
        exit(1);
        return init_dde_->Tbase();  // to make compiler happy
    }
}
const ChebyCoeff& DDC::Sbase() const {
    if (main_dde_)
        return main_dde_->Sbase();
    else if (init_dde_)
        return init_dde_->Sbase();
    else {
        std::cerr << "Error in DDC::Sbase(): Sbase is currently undefined" << std::endl;
        exit(1);
        return init_dde_->Sbase();  // to make compiler happy
    }
}

std::shared_ptr<DDE> DDC::newDDE(const std::vector<FlowField>& fields, const DDCFlags& flags) {
    std::shared_ptr<DDE> dde(new DDE(fields, flags));
    return dde;
}

std::shared_ptr<DDE> DDC::newDDE(const std::vector<FlowField>& fields, const std::vector<ChebyCoeff>& base,
                                 const DDCFlags& flags) {
    std::shared_ptr<DDE> dde(new DDE(fields, base, flags));
    return dde;
}

}  // namespace chflow