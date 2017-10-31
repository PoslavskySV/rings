//g++ -g -O2 factor_zp_ntl_flint.cpp -o factor_zp_ntl_flint -lntl -lgmp -lflint -lm

#include <NTL/ZZ_pXFactoring.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mod_poly.h>
#include <cstdio>

NTL_CLIENT

struct FlintZZ {
   fmpz_t value;

   explicit
   FlintZZ(const ZZ& a)
   {
      long n = NumBits(a);
   
      fmpz_init(value);
   
      for (long i = 0; i < n; i++)
         if (bit(a, i)) fmpz_setbit(value, i);
   
      if (a < 0)
         fmpz_neg(value, value);

   }

   ~FlintZZ() { fmpz_clear(value); }

};

struct FlintZZ_pX {
   fmpz_mod_poly_t value;

   explicit
   FlintZZ_pX(const ZZ_pX& a)
   {
      long da = deg(a);
      FlintZZ f_p(ZZ_p::modulus());
      fmpz_mod_poly_init2(value, f_p.value, da+1);

      for (long i = 0; i <= da; i++) {
         FlintZZ f_c(rep(a[i]));
         fmpz_mod_poly_set_coeff_fmpz(value, i, f_c.value);
      }
   }

   ~FlintZZ_pX() { fmpz_mod_poly_clear(value); }
};


struct FlintFac {
   fmpz_mod_poly_factor_t value;

   FlintFac() { fmpz_mod_poly_factor_init(value); }
   ~FlintFac() { fmpz_mod_poly_factor_clear(value); }

};

int main()
{
   setbuf(stdout, NULL);

   cout << "degree" << "\t" << "NTLTime" << "\t" << "FlintTime" << "\n";

   for (long l = 4; l <= 12; l += 1) {
      for (long idx = 0; idx < 10; idx++) {
         long degree = (1 << l) + idx * (((1 << (l + 1)) - (1 << l)) / 10);
         
         SetSeed((ZZ(l) << 64) + ZZ(degree));

         ZZ p(17);
         ZZ_p::init(p);

         ZZ_pX ntl_polynomial;
         ntl_polynomial.SetLength(degree + 1);
         for(int iCoef = 0; iCoef <= degree; ++iCoef)
            ntl_polynomial[iCoef] = iCoef;
         ntl_polynomial[0] = 1;
         ntl_polynomial.normalize();
         MakeMonic(ntl_polynomial);

         Vec< Pair<ZZ_pX, long> > ntl_factors;


         static int nIterations = 3;
         double t;

         t = GetTime();
         for (int i = 0; i < nIterations; i++)
            CanZass(ntl_factors, ntl_polynomial);
         double NTLTime = GetTime() - t;
   
         FlintZZ_pX flint_polynomial(ntl_polynomial);
         FlintFac flint_factors;

         t = GetTime();
         for (int i = 0; i < nIterations; i++)
            fmpz_mod_poly_factor_kaltofen_shoup(flint_factors.value, flint_polynomial.value);
         double FlintTime = GetTime()-t;

         cout << degree << "\t" << NTLTime << "\t" << FlintTime << "\n";
      }
   }
}


