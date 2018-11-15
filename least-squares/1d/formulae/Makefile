TARGETS := $(patsubst %.tex,%.png,$(wildcard *.tex))

all: $(TARGETS)
#	sudo cp *png /mnt/tmp/habr/

$(TARGETS): %.png: %.tex
	rubber $<
	dvisvgm --no-fonts --scale=1.5 `basename $< tex`dvi
	convert -density 85 `basename $< tex`svg $@
	rubber --clean $<

