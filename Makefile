
## comment this out if you need a different version of R, 
## and set set R_HOME accordingly as an environment variable
R_HOME := 		$(shell R RHOME)

## include headers and libraries for R 
RCPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --cppflags) -I/usr/include/gsl -lgsl -lgslcblas -lm
RLDFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --ldflags) -L/usr/lib -lgsl -lgslcblas -lm

## include headers and libraries for Rcpp interface classes
RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
ARMAINCL := 		$(shell echo 'RcppArmadillo:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)
GSLLIB   :=             $(shell gsl-config --libs)
GSLINC   :=             $(shell gsl-config --cflags)

## OpenMP
OPENMPFLAGS :=		-fopenmp

c_sources := 		$(wildcard *.c)


cpp_sources := 		$(wildcard *.cpp)


model.so :		$(c_sources) $(cpp_sources)
			PKG_CPPFLAGS="$(RCPPFLAGS) $(GSLLIB) $(GSLINC) $(RCPPINCL) $(ARMAINCL) $(OPENMPFLAGS)" PKG_LIBS="$(RLDFLAGS) $(RCPPLIBS) $(OPENMPFLAGS)" R CMD SHLIB $^

clean : 		
			/bin/rm *.o *.so 
