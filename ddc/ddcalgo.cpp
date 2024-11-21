/**
 *
 * Original author: Duc Nguyen
 */

#include "modules/ddc/ddcalgo.h"

namespace chflow {
    
void DDCAlgo::advance(vector<FlowField>& fieldsn, int Nsteps) {
    // This calculation follows Peyret section 4.5.1(b) pg 131.
    const int J = order_ - 1;
    vector<FlowField> rhs(nse_->createRHS(fieldsn));  // Number of fields and number of RHS's can be different
    int len = rhs.size();
    fields_[0] = fieldsn;
    // Start of time stepping loop

    for (int step = 0; step < Nsteps; ++step) {
        // Calculate nonlinearity, includes dealiasing if applicable
        if (order_ > 0) {
            nse_->nonlinear(fields_[0], nonlf_[0]);
        }

        // Add up multistepping terms of linear and nonlinear terms
        for (int l = 0; l < len; ++l) {
            rhs[l].setToZero();  // RHS must be zero before sum over multistep loop
            for (int j = 0; j < order_; ++j) {
                const Real a = -alpha_[j] / flags_.dt;
                const Real b = -beta_[j];
                rhs[l].add(a, fields_[j][l], b, nonlf_[j][l]);
            }
        }
        // Solve the implicit problem
        nse_->solve(fields_[J], rhs);
        // The solution is currently stored in fields_[J]. Shift entire fields and nonlf vectors
        // to move it into fields_[0]. Ie shift fields_[J] <- fields_[J-1] <- ... <- fields_[0] <- fields_[J]
        for (int j = order_ - 1; j > 0; --j) {
            for (int l = 0; l < numfields_; ++l) {
                swap(nonlf_[j][l], nonlf_[j - 1][l]);
                swap(fields_[j][l], fields_[j - 1][l]);
            }
        }
        #ifdef FREESLIP
        cout <<"Step "<<step<<"/"<<Nsteps<<endl;
        // add boundary conditions
        freeslipBC(fields_[0]);
        #endif

        t_ += flags_.dt;

        // printStack();
        //*flags_.logstream << "} Multistep::advance(...) step == " << step << " }" <<endl;
        if (nse_->taskid() == 0) {
            if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll)
                *flags_.logstream << t_ << ' ' << flush;
            else if (flags_.verbosity == PrintTicks)
                *flags_.logstream << '.' << flush;
        }
    }  // End of time stepping loop

    fieldsn = fields_[0];  // update velocity

    if (nse_->taskid() == 0)
        if (flags_.verbosity == PrintTime || flags_.verbosity == PrintAll || flags_.verbosity == PrintTicks)
            *flags_.logstream << endl;
    return;
}

}  // namespace chflow