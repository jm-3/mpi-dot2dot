#include "algorithm.h"
#include <omp.h>

int start_TRs_search (Dot_input* param) {
	int window_length, window_length_max, window_index, window_end, max_length_flag, 
		jumps_limit=0, gaps_limit=0, biggest_full_length=0;
	//MATCH_ARRAY_TYPE minThresholdByPercentage, minThresholdByGaps, minThreshold;
	float minThresholdByPercentage, minThresholdByGaps, minThreshold;	
	struct dot_matrix *m = param->matrix;
	struct config *cfg = param->config_params;
	result_findTR *current_result;	
	TRs_Result_Bundle *trs_current_bundle, *last_tandem_found, *previous_window_tandem;
	int error=0;

	if ((cfg->fvalue < 0) || (cfg->fvalue > 1)) {
		printf ("MinMatch ranges in [0,1]\n");
		return 1;
	}
	if (cfg->gvalue < 0) { 
		printf ("max_gaps cannot be negative\n");
		return 1; 
	} 
	else {
		gaps_limit=cfg->gvalue;
	}

	if (cfg->xvalue < 1) { /* max_length not inserted in command line */
		max_length_flag = 0;
	} 
	else { /* max_length inserted */
		if (cfg->xvalue < cfg->nvalue) { /* if lower than min_length */
			printf ("max_length < min_length not allowed\n");
			return 1;
		}
		max_length_flag = 1;
	}
	if (cfg->nvalue < 1) { 
		printf ("motif length cannot be lower than 1\n");
		return 1;
	} 
	else { 
		window_length = cfg->nvalue; 
	}

	window_index=0;   
	window_end=m->sequence_len; 


	if (cfg->jvalue < 0) {
		printf ("MaxInsert cannot be negative\n");
		return 1;
	}

	if ((current_result = init_result_struct()) == NULL) {
		perror ("Error in start_TRs_search() for current_result initialization\n");
		return 1;
	}

	if (( trs_current_bundle = init_TRs_Bundle(RESIZE_TR_BUNDLE_AMOUNT, RESIZE_TR_MOTIFS_AMOUNT) ) == NULL) {
			perror("Error in initializing trs_current_bundle\n");
			return 1;
	}

	if (( previous_window_tandem = init_TRs_Bundle( 1, RESIZE_TR_MOTIFS_AMOUNT ) ) == NULL) {
			perror("Error in initializing previous_window_tandem\n");
			return 1;
	}

	if (( last_tandem_found = init_TRs_Bundle( 1, RESIZE_TR_MOTIFS_AMOUNT ) ) == NULL) {
			perror("Error in initializing trs_current_bundle\n");
			return 1;
	}

	#pragma omp parallel for reduction(+:error) shared(window_end)
	for(window_index = 0; window_index < window_end; window_index++ ) {
#ifdef DEBUG_ALG
		printf("TRS WINDOW INDEX STARTED:%d\n", window_index);
#endif
		if (max_length_flag && ( cfg->xvalue <= ((m->sequence_len - window_index)/2) )) {
			window_length_max = cfg->xvalue;			
		} else {
			window_length_max = (m->sequence_len - window_index)/2;			
		}		
		
		
		biggest_full_length = 0;
		while( window_length <= window_length_max ) {			
#ifdef DEBUG_ALG
			printf("\tTRS WINDOW LENGTH STARTED:%d\n", window_length);
#endif					
			//minThresholdByPercentage = (MATCH_ARRAY_TYPE) cfg->fvalue * window_length;
			minThresholdByPercentage = (float) cfg->fvalue * window_length;
			//minThresholdByGaps = (gaps_limit >= window_length) ? ((MATCH_ARRAY_TYPE) window_length-1) : ((MATCH_ARRAY_TYPE) window_length - gaps_limit);
			minThresholdByGaps = (gaps_limit >= window_length) ? ((float) window_length-1) : ((float) window_length - gaps_limit);					
					
			if ((cfg->fvalue == 1) && (gaps_limit > 0)) { /* only gaps parameter */
				minThreshold = minThresholdByGaps;
			} else if ((gaps_limit == 0) && (cfg->fvalue < 1)) { /* only percentage parameter */
				minThreshold = minThresholdByPercentage;
			} else { /* both or no parameter: percentage and gaps */
				minThreshold = MAX(minThresholdByPercentage, minThresholdByGaps); 
			}

#ifdef DEBUG_ALG
			printf("\tMINTHRESHOLD %f\n", minThreshold);
#endif	
			/* Jump limit is at most window_length or cfg->jvalue */
			if (cfg->jvalue > window_length) { jumps_limit = window_length; } else { jumps_limit = cfg->jvalue; }

			if (findTandemRepeats(window_length, window_index, m, minThreshold, jumps_limit, current_result)) {
				printf ("Error in start_TRs_search() for current_result\n");
				error = 1;
				window_end = 0; /* terminate all threads */
			}
			switch (current_result->result_code) {
				case (NO_TR_FOUND) : {
					window_length++;
					break;
				}
				case (TR_FOUND) : { 
					window_length ++;
					/* Put all the TRs with the same window_index in the TRs_Bundle_list and check which is the highest length */
					if (current_result->resulted_TR->TRs_found[0].full_length > biggest_full_length) biggest_full_length = current_result->resulted_TR->TRs_found[0].full_length;
					if (insert_TRresult_inBundle(trs_current_bundle, current_result->resulted_TR, RESIZE_TR_BUNDLE_AMOUNT, RESIZE_TR_MOTIFS_AMOUNT)) {
						perror("Error in inserting the current result in the TRs bundle\n");
						error = 1;
					}
#ifdef DEBUG_ALG
					printf("\tTRFOUND at %d INSERTED\n", window_index);
#endif
					break;
				}
				default : {
					printf ("Unknown case!\n");
					error = 1;
					window_end = 0; /* terminate all threads */
				}
			}
			/* Reset current_result  */
			reset_result_struct(current_result);
			
		}
		
		/* Apply Expansion Filter to the TRs list of the same zone (same Window_index)*/
		/*last_tandem_found is NULL or has a valid TR to compare to the previous one */
		expansion_filter(trs_current_bundle, last_tandem_found, biggest_full_length, cfg->tollerance);
#ifdef DEBUG_ALG
		printf("\tTRS_BUNDLE of %d FILTERED\n", window_index);
#endif
		
		if ( last_tandem_found->motif_lengths_offset != 0 ) {/* Tandem found */
		  /*Insert the tandem as it is the first found	 */
			if (previous_window_tandem->motif_lengths_offset == 0) {
#ifdef DEBUG_ALG
				printf("\tTRS WINDOW INDEX PREVIUOS WITH 0 COPY NUMBER:%d\n", window_index);
#endif
				if (insert_TRresult_inBundle(param->TRs_bundle , last_tandem_found, RESIZE_TRS_AMOUNT, RESIZE_MOTIFS_AMOUNT)) {
					perror("Error in inserting the current result in the TRs bundle\n");
					error = 1;
					window_end = 0; /* terminate all threads */
				};
				if ( copy_TRs_Bundle(last_tandem_found, previous_window_tandem) ) {
					perror("Error in copying last_tandem_found to previous_window_tandem in start_TRs_search\n");
					error = 1;
					window_end = 0; /* terminate all threads */
				};
				//window_index++;
			} else {
				
				if (isLastIncluded(&(previous_window_tandem->TRs_found[0]), &(last_tandem_found->TRs_found[0]))) {
					/* included Tandems */
#ifdef DEBUG_ALG
					printf("\tTRS WINDOW INDEX INCLUDED:%d\n", window_index);
#endif
					if (last_tandem_found->TRs_found[0].purity_percentage > previous_window_tandem->TRs_found[0].purity_percentage/*(last_tandem_found->period > previous_window_tandem->period)*/) { /* if false it has found an included TR which is the same but fragmented */
						switch (m->mask[window_index]) {
							case (UNCHECKED) : {
								/* tandem has not been checked before */ 				
								if (insert_TRresult_inBundle(param->TRs_bundle , last_tandem_found, RESIZE_TRS_AMOUNT, RESIZE_MOTIFS_AMOUNT)) {
									perror("Error in inserting the current result in the TRs bundle\n");
									error = 1;
									window_end = 0; /* terminate all threads */
								};
								//window_index++;
								break;											
							}
							case (CHECKED) : {
								/* it is a part of an other tandem found before */
								//window_index++;
								break;
							}
							default : { break; }
						}
					}
				} else { /* intersected Tandems */		
#ifdef DEBUG_ALG
					printf("\tTRS WINDOW INDEX INTERSECTED:%d\n", window_index);				
#endif
					if (insert_TRresult_inBundle(param->TRs_bundle , last_tandem_found, RESIZE_TRS_AMOUNT, RESIZE_MOTIFS_AMOUNT)) {
						perror("Error in inserting the current result in the TRs bundle\n");
						error = 1;
						window_end = 0; /* terminate all threads */
					};
					/*reset_TRs_Bundle(previous_window_tandem);*/
					if ( copy_TRs_Bundle(last_tandem_found, previous_window_tandem) ) {
						perror("Error in copying last_tandem_found to previous_window_tandem in start_TRs_search\n");
						error = 1;
						window_end = 0; /* terminate all threads */
					}
					//window_index++;
				}
			}
		}
		reset_TRs_Bundle(last_tandem_found);
		reset_TRs_Bundle(trs_current_bundle);
		window_length=cfg->nvalue;
	
	}
	destroy_TRs_Bundle(&last_tandem_found);
	destroy_TRs_Bundle(&previous_window_tandem);
	destroy_result_struct(&current_result);
	destroy_TRs_Bundle(&trs_current_bundle);
	return error;
}