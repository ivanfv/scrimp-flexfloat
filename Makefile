all:
	gcc -O3 ./scrimp_ff/scrimp_ff.c -o ./scrimp_ff/scrimp_ff -L./flexfloat/build -lflexfloat -lm -fopenmp
clean:
	rm ./scrimp_ff/scrimp_ff
random_similarity:
	./scrimp_ff/scrimp_ff random_similarity.txt 50 4 1
random_anomaly:
	./scrimp_ff/scrimp_ff random_anomaly.txt 50 4 1
