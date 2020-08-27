#include <stdio.h>
#include <pthread.h>
#include <unistd.h>
#include <stdlib.h>

#define NUMTHREADS 20

pthread_mutex_t mutex;
pthread_cond_t cond[NUMTHREADS];
int contador=0;


void *funcion(void * params){
	int id = *((int *) params);

	srand(id);
	printf("inicia hebra %d\n", id);
while(1){
	pthread_mutex_lock(&mutex);
	printf("hebra %d adquiere lock\n",id);
	/*while(id != contador){
		printf("hebra %d se bloquea esperando cond\n",id);
		pthread_cond_wait(&cond, &mutex);
		printf("hebra %d liberada por cond\n",id);
	}*/
	if(id != contador){
		printf("hebra %d se bloquea esperando cond\n",id);
		pthread_cond_wait(&cond[id], &mutex);
		printf("hebra %d liberada por cond\n",id);
	}
	printf("hebra %d lee\n",id);
	sleep((rand()%10)+1);
	contador = (contador + 1)%NUMTHREADS;
	//pthread_cond_broadcast(&cond);
	pthread_cond_signal(&cond[(id+1)%NUMTHREADS]);
	printf("hebra %d libera lock\n",id);
	pthread_mutex_unlock(&mutex);
}
}



int main(){
	pthread_t hebras[NUMTHREADS];
	int i, ids[NUMTHREADS];

	pthread_mutex_init(&mutex, NULL);
	for(i=0; i<NUMTHREADS; i++)
		pthread_cond_init(&cond[i], NULL);

	for(i=0; i<NUMTHREADS; i++){
		ids[i] = i;
		pthread_create(&hebras[i],NULL, &funcion, &ids[i]);
	}

	for(i=0; i<NUMTHREADS; i++){
		pthread_join(hebras[i], NULL);
	}

	return 0;
}
