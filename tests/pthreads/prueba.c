#include <stdio.h>
#include <pthread.h>
#include <unistd.h>

pthread_mutex_t mutex;
pthread_cond_t cond;
int contador=0;

int NUMTHREADS=10;

void *funcion(void * params){
	int id = *((int *) params);

	printf("inicia hebra %d\n", id);
while(1){
	pthread_mutex_lock(&mutex);
	printf("hebra %d adquiere lock\n",id);
	while(id != contador){
		printf("hebra %d se bloquea esperando cond\n",id);
		pthread_cond_wait(&cond, &mutex);
		printf("hebra %d liberada por cond\n",id);
	}
	sleep(2);
	contador = (contador + 1)%NUMTHREADS;
	pthread_cond_broadcast(&cond);
	//pthread_cond_signal(&cond);
	printf("hebra %d libera lock\n",id);
	pthread_mutex_unlock(&mutex);
}
}



int main(){
	pthread_t hebras[NUMTHREADS];
	int i, ids[NUMTHREADS];

	pthread_mutex_init(&mutex, NULL);
	pthread_cond_init(&cond, NULL);

	for(i=0; i<NUMTHREADS; i++){
		ids[i] = i;
		pthread_create(&hebras[i],NULL, &funcion, &ids[i]);
	}

	for(i=0; i<NUMTHREADS; i++){
		pthread_join(hebras[i], NULL);
	}

	return 0;
}
