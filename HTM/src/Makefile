#Makefile for building and running the LSMR model on a Linux or Mac OSX
EXEC=HTM
TPL=$(EXEC).tpl
DAT=$(EXEC).dat
ARG=


all: $(EXEC) $(EXEC).rep dust

$(EXEC): $(TPL)
	admb $(EXEC)

$(EXEC).rep: $(TPL) $(DAT)
	./$(EXEC) -ind $(DAT) $(ARG)

mcmc: runmc mceval dust

runmc: $(EXEC) $(EXEC).rep
	./$(EXEC) -ind $(DAT) $(ARG) -mcmc 500000 -mcsave 200 -nosdmcmc

mceval: $(EXEC).psv
	./$(EXEC) -ind $(DAT) $(ARG) -mceval



dust:
	rm -f *.log *.rpt *.htp admodel.* variance *.bar *.mcm *.[bpr][0123456789]*

clean:
	rm -f $(EXEC) *.eva *.log *.rpt *.htp *.cor *.par *.r* *.p0* *.b*
	rm -f *.rep *.bar *.psv *.std $(EXEC).cpp admodel.* variance