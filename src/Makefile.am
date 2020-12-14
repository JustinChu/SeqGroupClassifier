bin_PROGRAMS = kc

kc_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

kc_CPPFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/vendor \
	-I$(top_srcdir)
		
kc_LDFLAGS = $(OPENMP_CXXFLAGS) \
	-lz \
	-pthread

kc_SOURCES = KmerClassifier.cpp