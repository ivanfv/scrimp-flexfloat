all:
	gcc -O3 ./scrimp/scrimp.c -o ./scrimp/scrimp -L./flexfloat/build -lflexfloat -lm -fopenmp
clean:
	rm ./scrimp/scrimp
random_similarity:
	./scrimp/scrimp random_similarity.txt 50 72 1
random_anomaly:
	./scrimp/scrimp random_anomaly.txt 50 72 1
