# Makefile for the maintainence of the GDAGsim library

NAME=gdagsim-03

backup:
	make rebuild
	cp $(NAME).tgz Backup/$(NAME).`date +"%Y-%m-%d"`.tgz

web:
	make rebuild
	scp $(NAME).tgz Doc/gdag.pdf @finan:public_html/software/gdagsim/

rebuild:
	rm -f $(NAME).tgz
	make $(NAME).tgz


$(NAME).tgz:
	cd Doc ; make
	make clean
	tar cfz $(NAME).tgz $(NAME)

clean:
	rm -f *~ *.tgz
	cd $(NAME) ; rm -f *~
	cd $(NAME)/src ; make clean
	cd $(NAME)/examples ; make clean

full-clean:
	make clean
	cd Doc ; make full-clean
	cd Examples ; make clean
 
work:
	make full-clean
	scp -r Makefile Doc Examples $(NAME) @bayes.ncl.ac.uk:src/gdag

home:
	make full-clean
	scp -r Makefile Doc Examples $(NAME) @${REMOTEHOST}:src/gdag




# eof

