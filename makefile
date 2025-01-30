Main = ./projekt8.c
Sources = MyNumerics/my_numerics.c
Output = projekt8

Libs = -lm -lgsl -lgslcblas 

Flags = -Wall -Wextra

#All-Befehl
all: compile run

compile:
	gcc $(Main) $(Sources) -o $(Output) $(Libs) $(Flags)

compile_debug:
	gcc $(Main) $(Sources) -o $(Output) -g -fsanitize=address $(Libs) $(Flags)


run:
	./$(Output)

run_debug: compile_debug
	./$(Output)

clean:
	rm -f $(Output)
	rm -f *.txt
	rm -f *.csv

plot:
	python3 plots.py