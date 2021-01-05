bin_PROGRAMS = seqgc

seqgc_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

seqgc_CPPFLAGS = -I$(top_srcdir)/src \
	-I$(top_srcdir)/vendor \
	-I$(top_srcdir)
		
seqgc_LDFLAGS = $(OPENMP_CXXFLAGS) \
	-lz \
	-pthread

seqgc_SOURCES = SeqGC.cpp