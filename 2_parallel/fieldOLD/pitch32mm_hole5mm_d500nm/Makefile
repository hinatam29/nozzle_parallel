# ----------------------------------------------------
#     MAKEFILE FOR THE CALUCURATION OF FLOW FIELD
# ----------------------------------------------------
#
CAL = main.o geo.o initial.o bound.o advect.o cipcsl3x.o cipcsl3y.o TECSIAx.o TECSIAy.o result.o pfield_sor.o pfield_bicg.o proj.o proj_siax.o proj_siay.o TECVIAx.o TECVIAy.o crrctx.o crrcty.o diffusion.o advVOF.o grav.o prop.o dfunc.o csf.o average.o average2.o nonadv.o kzoku.o kinput.o potential.o ej.o charge.o eforce.o fixvf0.o elec.o difrho.o heviside.o viscous.o pserch.o bound1.o psearch.o intraction.o rancal.o pcount.o diffusion_liq.o diffusion_gas.o diffusion1.o

#
tpfexp :  ${CAL}
	gfortran -fopenmp -o go.out ${CAL}
#
.f.o :
#
#	ifort -c -O3 -autodouble -fpe0 -traceback -check all -g $*.f
#	ifort -c -O3 -autodouble $*.f
	gfortran -c -O3 -fopenmp -fdefault-real-8 -std=legacy -fallow-argument-mismatch $*.f
#
${CAL}: subcom.inc
clean: 
