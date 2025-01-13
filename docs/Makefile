all:  clean build open

bo: build open

build: $(wildcard _build/src/*)
	jupyter-book build .

clean:
	rm -rf _build

open:
	open _build/html/index.html 