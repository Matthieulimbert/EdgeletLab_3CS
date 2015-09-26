include ../../Makefile.include

TARGETS = 
SUBDIRS = Networks

all : $(TARGETS)
	for d in $(SUBDIRS); do \
		(if [ -d $$d ]; then cd $$d && ${MAKE} all; fi); \
	done

clean :
	for d in $(SUBDIRS); do \
		(if [ -d $$d ]; then cd $$d && ${MAKE} clean; fi); \
	done; \
	rm -f $(TARGETS)
