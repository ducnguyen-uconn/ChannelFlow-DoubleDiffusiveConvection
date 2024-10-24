/**
 * Main interface to handle DNS of DDC
 *
 * Original author: Duc Nguyen
 */

#ifndef DDC_H
#define DDC_H

#include "channelflow/dns.h"
#include "channelflow/dnsalgo.h"
#include "modules/ddc/ddcflags.h"
#include "modules/ddc/dde.h"
#include "modules/ddc/ddcdsi.h"

namespace chflow {

int field2vector_size(const FlowField& u, const FlowField& temp, const FlowField& salt);

/** \brief Turn the three flowfields for velocity and temperature into one Eigen vector
 * \param[in] u velocity field
 * \param[in] temp temperature field
 * \param[in] salt salinity field
 * \param[in] x vector for the linear algebra
 *
 * The vectorization of u is analog to the field2vector, temparture and salinity are piped entirely
 * into the vector (a single independent dimensions)
 */
void field2vector(const FlowField& u, const FlowField& temp, const FlowField& salt, Eigen::VectorXd& x);

/** \brief Turn  one Eigen vector into the three flowfields for velocity, temperature and salinity
 * \param[in] x vector for the linear algebra
 * \param[in] u velocity field
 * \param[in] temp temperature field
 * \param[in] salt salinity field
 *
 * The vectorization of u is analog to the field2vector, temperature and salinity are piped entirely
 * into the vector (a single independent dimension)
 */
void vector2field(const Eigen::VectorXd& x, FlowField& u, FlowField& temp, FlowField& salt);

/** \brief extension of the DNSFlags class
 *
 * VEDNSFlags class, holds all additional parameters for viscoelastic fluids
 */

/** \brief wrapper class of DNSAlgorithm and DDC
 *
 *
 */
class DDC : public DNS {
   public:
    //     DDC ();
    //     DDC (const DDC & ddc);
    DDC(const std::vector<FlowField>& fields, const DDCFlags& flags);
    //     DDC (const vector<FlowField> & fields, const vector<ChebyCoeff> & base,
    // 	    const DDCFlags & flags);

    virtual ~DDC();

    DDC& operator=(const DDC& ddc);

    //     virtual void advance (vector<FlowField> & fields, int nSteps = 1);
    //
    //     virtual void reset_dt (Real dt);
    //     virtual void printStack () const;

    const ChebyCoeff& Ubase() const;
    const ChebyCoeff& Wbase() const;
    const ChebyCoeff& Tbase() const;
    const ChebyCoeff& Sbase() const;

   protected:
    std::shared_ptr<DDE> main_dde_;
    std::shared_ptr<DDE> init_dde_;

    std::shared_ptr<DDE> newDDE(const std::vector<FlowField>& fields, const DDCFlags& flags);
    std::shared_ptr<DDE> newDDE(const std::vector<FlowField>& fields, const std::vector<ChebyCoeff>& base,
                                const DDCFlags& flags);
};


}  // namespace chflow
#endif  // DDC_H