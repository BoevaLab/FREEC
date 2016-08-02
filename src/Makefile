#
# FREEC project Makefile wrapper
#
# Institut Curie, France
#
# Eric Viara, Valentina Boeva January 2013
#

MAKEFILE = Makefile.freec

all:
%:
	@rm -f depend.mk
	@touch depend.mk
	$(MAKE) -f $(MAKEFILE) init
	$(MAKE) -f $(MAKEFILE) depend
	$(MAKE) -f $(MAKEFILE) $@

clean:
	@rm -f depend.mk
	@touch depend.mk
	$(MAKE) -f $(MAKEFILE) $@

