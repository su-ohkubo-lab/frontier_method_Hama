#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_ID 10000000
#define MAX_LEN_STR 4096

#define TRUE 1
#define FALSE 0
#define IS_CUTOFF_USED TRUE
#define SEC2NSEC (1 * 1000 * 1000 * 1000)
//#define IS_CUTOFF_USED FALSE

typedef long long int int64; // 64 bit integer

// Declare global variables from the input files.
// *** NOTE:
// *** ori_event: The original events in the input file
// *** event: Events with the same state change (but different rate constants)
static int NUM_ORI_EVENT; // The number of ori_events
static int NUM_EVENT; // The number of events (without state-change dupilication)
static int DIM_STATE; // The dimension of state (the number of components)
static int NUM_CONSTANT; // The number of constants in the input file
//--- The followings are used for the correspondence between ori_events and events.
static int **ORI_EVENT;
static int **EVENT;
static double *ORI_RATE;
static int **ORI_INDEX4RATE;
static int **ORI_INDEX4CONST;
static int *LIST_LEN_E2ORI;
static int *LIST_START_INDEX_E2ORI;
//--- The followings are needed for the time-evolution process.
static int ID_DIAGONAL_EVENT;
static double *PRUNING_COEFF_POS;
static double *PRUNING_COEFF_NEG;
static double FINAL_TIME;
static int *INITIAL_STATE;
static int *TARGET_STATE;
static double *CONSTS;

// Declare working global variables.
static double *ORI_COEFF;
static int **ORI_STATE;
static double *ORI_EVENT_CONST_RATE;
static double *COEFF;
static int **STATE;
// --- Only for check
int64 MAXIMUM_ID_NUMBER;

// Declare subroutines for the settings.
void load_event_file(char *filename);
void load_parameter_file(char *filename);
void initialize_settings_for_dual_comp(void);
void set_pruning_coeffs(void);
void set_state_from_str(int *state, char *str);
int is_pruned(int *current, int *target, int length);
int is_out_of_range(int *current);

// ---- The followings are needed for Crank-Nicolson method
static int **LIST_EVENTS_IN_DIAGONAL_RELAY_PATH;
static int LENGTH_OF_DIAGONAL_RELAY_PATH;
void make_2nd_order_event(void);
void free_memories_allocated_in_make_2nd_order_event(void);

// Declare subroutines for the time-evolution.
double evaluate_with_Euler(double dt, int max_M);
double evaluate_with_Heun(double dt, int max_M);
double evaluate_with_resolvent(double dt, int max_M);
double evaluate_with_Crank_Nicolson_simple(double dt, int max_M);
double evaluate_with_Crank_Nicolson(double dt, int max_M);
double calculate_transition_matrix_element(int event, int *state);
int is_valid_next_state(int event, int *state);
double occur_event_Euler(int event, int *state, double dt_M);
double occur_event_Heun_2nd_process(int event, int *state, double dt_M);
double occur_event_resolvent(int event, int *state, double dt_M);
double occur_event_resolvent_2nd_order(int event, int *state, double dt_M);
double evaluate_diagonal_factor(int *state, double dt_M);

// Declare subroutines for file open and memory allocations.
FILE *f_open(const char *path, const char *mode);
void *clear_alloc(size_t size, size_t nobj);
void **clear_alloc_2d( size_t size, size_t row, size_t column);
#define free_2d(p) {free(p[0]);free(p);}
void free_memories_allocated_in_load_event_file(void);

void diff_timespec(struct timespec *result, struct timespec *end,
                   struct timespec *begin)
{
  result->tv_sec  = end->tv_sec  - begin->tv_sec;
  result->tv_nsec = end->tv_nsec - begin->tv_nsec;
 
  if (result->tv_nsec < 0) {
    result->tv_sec--;
    result->tv_nsec += SEC2NSEC;
  }
}

int main(int argc, char *argv[]){
  if(argc != 3){
    fprintf(stderr, "Usage: ./app_dual_computation event_input_file parameter_input_file\n" );
    fprintf(stderr, "example: ./app_dual_computation tmp_event.txt tmp_parameters.txt\n" );
    exit(1);
  }
  struct timespec tv, tv2, result;

  // Load the input files (event and parameter).
  load_event_file(argv[1]);
  INITIAL_STATE = (int*)clear_alloc(sizeof(int), DIM_STATE);
  TARGET_STATE = (int*)clear_alloc(sizeof(int), DIM_STATE);
  CONSTS = (double*)clear_alloc(sizeof(double), NUM_CONSTANT);
  load_parameter_file(argv[2]);

  // Initialize settings.
  initialize_settings_for_dual_comp();
  make_2nd_order_event();
  
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tv);
  // Perform the main calculation.  // 90 45 30  // 24  12  8
  double comp = evaluate_with_Euler(FINAL_TIME, 5);
  //double comp = evaluate_with_Heun(FINAL_TIME, 45);
  // double comp = evaluate_with_resolvent(FINAL_TIME, 20);
  //double comp = evaluate_with_Crank_Nicolson(FINAL_TIME, 30);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tv2);
  diff_timespec(&result, &tv2, &tv);
  printf("Execution Time = %lf [s]\n", result.tv_sec + (double) result.tv_nsec / SEC2NSEC);
  printf("%.15e\n", comp);

  // --- Only for check
  //printf("MAXIMUM_ID_NUMBER: %lld\n", MAXIMUM_ID_NUMBER);
  
  // Free the allocated memories.
  free_memories_allocated_in_make_2nd_order_event(); 
  free(ORI_EVENT_CONST_RATE);
  free_2d(ORI_STATE);
  free_2d(STATE);
  free(ORI_COEFF);
  free(COEFF);
  free_memories_allocated_in_load_event_file();
  free(PRUNING_COEFF_NEG);
  free(PRUNING_COEFF_POS);
  free(CONSTS);
  free(TARGET_STATE);
  free(INITIAL_STATE);
  
  return 0;  
}

// -----------------------------------------------------------
// Subroutines 
// -----------------------------------------------------------
void load_event_file(char *filename){
  FILE *fp;
  char buf[MAX_LEN_STR];
  char str_state_changes[MAX_LEN_STR];
  char str_state_changes_previous[MAX_LEN_STR];
  char str_rate_coeff[MAX_LEN_STR];
  char str_rate_factors[MAX_LEN_STR];
  char str_rate_constants[MAX_LEN_STR];
  fp = f_open(filename, "r");
  // -- Skip the first line.
  fgets(buf, sizeof(buf), fp);
  // -- Read parameters.
  fscanf(fp, "%d %d %d %d\n", &NUM_ORI_EVENT, &NUM_EVENT, &DIM_STATE, &NUM_CONSTANT);
  // -- Skip the third line.
  fgets(buf, sizeof(buf), fp);
  // Allocate memories.
  ORI_EVENT = (int**)clear_alloc_2d(sizeof(int), NUM_ORI_EVENT, DIM_STATE);
  EVENT = (int**)clear_alloc_2d(sizeof(int), NUM_EVENT, DIM_STATE);
  ORI_RATE = (double*)clear_alloc(sizeof(double), NUM_ORI_EVENT);
  ORI_INDEX4RATE = (int**)clear_alloc_2d(sizeof(int), NUM_ORI_EVENT, DIM_STATE);
  ORI_INDEX4CONST = (int**)clear_alloc_2d(sizeof(int), NUM_ORI_EVENT, NUM_CONSTANT);
  LIST_LEN_E2ORI = (int*)clear_alloc(sizeof(int), NUM_EVENT + 1); // '+1' is needed for the case in which the diagonal part is zero.
  LIST_START_INDEX_E2ORI = (int*)clear_alloc(sizeof(int), NUM_EVENT);
  for(int i = 0; i < NUM_EVENT+1; i++){
    LIST_LEN_E2ORI[i] = 0;
  }
  int e = 0;
  int ori_e = 0;
  int is_there_a_diagonal_event = FALSE;
  strcpy(str_state_changes_previous, "");
  while(fgets(buf, sizeof(buf), fp) != NULL){
    char str_tmp[4096];
    if(*buf == '\n' || buf[0] == ' ') continue;
    sscanf(buf,"%[^/] / %[^/] / %[^/] / %[^#]", str_state_changes, str_rate_coeff, str_rate_factors, str_rate_constants);
    {
      // Read state_change.
      if(strcmp(str_state_changes, str_state_changes_previous) != 0){
	// Set EVENT and ORI_EVENT.
	char tmp[4096];
	strcpy(tmp, str_state_changes);
	char *tp;
	int c = 0;
	int is_zero_change = TRUE;
	tp = strtok(tmp, " ,[]");
	ORI_EVENT[ori_e][c] = atoi(tp);
	EVENT[e][c] = atoi(tp);
	if(EVENT[e][c] != 0){
	  is_zero_change = FALSE;
	}
	c++;
	while(tp != NULL){
	  tp = strtok(NULL," ,[]");
	  if (tp != NULL){
	    ORI_EVENT[ori_e][c] = atoi(tp);
	    EVENT[e][c] = atoi(tp);
	    if(EVENT[e][c] != 0){
	      is_zero_change = FALSE;
	    }
	    c++;
	  }
	}
	// Set the correspondence from the event to the ori_event.
	LIST_START_INDEX_E2ORI[e] = ori_e;
	LIST_LEN_E2ORI[e]++;
	// Check whether there is an no-state-change event or not.
	if(is_zero_change == TRUE){
	  is_there_a_diagonal_event = TRUE;
	  ID_DIAGONAL_EVENT = e;
	}
	e++;
      }else{
	// Set only ORI_EVENT.
	char tmp[4096];
	strcpy(tmp, str_state_changes);
	char *tp;
	int c = 0;
	int is_zero_change = TRUE;
	tp = strtok(tmp, " ,[]");
	ORI_EVENT[ori_e][c] = atoi(tp);
	c++;
	while(tp != NULL){
	  tp = strtok(NULL," ,[]");
	  if (tp != NULL){
	    ORI_EVENT[ori_e][c] = atoi(tp);
	    c++;
	  }
	}
	// Set the correspondence from the event to the ori_event.
	LIST_LEN_E2ORI[e-1]++;
      }
      strcpy(str_state_changes_previous, str_state_changes);
      {
	// Read rate_coeff.
	ORI_RATE[ori_e] = atof(str_rate_coeff);
      }
      {
	// Read rate_factors.
	char *tp;
	int c = 0;
	tp = strtok(str_rate_factors, " ,[]");
	ORI_INDEX4RATE[ori_e][c] = atoi(tp);
	c++;
	while(tp != NULL){
	  tp = strtok(NULL," ,[]");
	  if (tp != NULL){
	    ORI_INDEX4RATE[ori_e][c] = atoi(tp);
	    c++;
	  }
	}
      }
      {
	// Read index for constant.
	char *tp;
	int c = 0;
	tp = strtok(str_rate_constants, " ,[]");
	ORI_INDEX4CONST[ori_e][c] = atoi(tp);
	c++;
	while(tp != NULL){
	  tp = strtok(NULL," ,[]");
	  if (tp != NULL){
	    ORI_INDEX4CONST[ori_e][c] = atoi(tp);
	    c++;
	  }
	}
      }
      ori_e++;
    }
  }
  fclose(fp);
  if(is_there_a_diagonal_event == FALSE){
    ID_DIAGONAL_EVENT = NUM_EVENT;
  }
}

void load_parameter_file(char *filename){
  FILE *fp;
  char buf[MAX_LEN_STR];
  fp = f_open(filename, "r");
  // -- Read the final time.
  fgets(buf, sizeof(buf), fp); // Skip a comment.
  fgets(buf, sizeof(buf), fp);
  sscanf(buf, "%lf ", &FINAL_TIME);
  // -- Read the initial state.
  fgets(buf, sizeof(buf), fp); // Skip a comment.
  fgets(buf, sizeof(buf), fp);
  set_state_from_str(INITIAL_STATE, buf);
  // -- Read the target state.
  fgets(buf, sizeof(buf), fp); // Skip a comment.
  fgets(buf, sizeof(buf), fp);
  set_state_from_str(TARGET_STATE, buf);
  // -- Read the parameters.
  fgets(buf, sizeof(buf), fp); // Skip a comment.
  for(int i=0; i<NUM_CONSTANT; i++){
    fgets(buf, sizeof(buf), fp);
    double comp;
    sscanf(buf, "%lf ", &comp);
    CONSTS[i] = comp;
  }
  fclose(fp);
}
  
void free_memories_allocated_in_load_event_file(void){
  free_2d(ORI_EVENT);
  free_2d(EVENT);
  free(ORI_RATE);
  free_2d(ORI_INDEX4RATE);
  free_2d(ORI_INDEX4CONST);
  free(LIST_LEN_E2ORI);
  free(LIST_START_INDEX_E2ORI);
}

void initialize_settings_for_dual_comp(void){
  /// Set coefficients for the pruning.
  set_pruning_coeffs();
  // Memory allocation.
  COEFF = (double*)clear_alloc(sizeof(double), MAX_ID);
  ORI_COEFF = (double*)clear_alloc(sizeof(double), MAX_ID);
  STATE = (int**)clear_alloc_2d(sizeof(int), MAX_ID, DIM_STATE);
  ORI_STATE = (int**)clear_alloc_2d(sizeof(int), MAX_ID, DIM_STATE);
  ORI_EVENT_CONST_RATE = (double*)clear_alloc(sizeof(double), NUM_ORI_EVENT);
}

void set_pruning_coeffs(void){
  // Find the maximum step in the positive and negative directions in each dimension.
  PRUNING_COEFF_POS = (double*)clear_alloc(sizeof(double), DIM_STATE);
  PRUNING_COEFF_NEG = (double*)clear_alloc(sizeof(double), DIM_STATE);
  for(int d = 0; d < DIM_STATE; d++){
    int max = 0;
    int min = 0;
    for(int e = 0; e < NUM_EVENT; e++){
      if(EVENT[e][d] > max){
	max = EVENT[e][d];
      }
      if(EVENT[e][d] < min){
	min = EVENT[e][d];
      }
    }
    if(max != 0){
      PRUNING_COEFF_POS[d] = 1.0/(double)max;
    }else{
      fprintf(stderr, "There is no event with positive direction in %d th dimension.\n", d);
    }
    if(min != 0){
      PRUNING_COEFF_NEG[d] = 1.0/(double)min;
    }else{
      fprintf(stderr, "There is no event with negative direction in %d th dimension.\n", d);
    }
  }
}

void set_state_from_str(int *state, char *str){
  char *tp;
  int d = 0;
  tp = strtok(str, " ,'\"");
  state[d] = atoi(tp);
  d++;
  while (tp != NULL){
    tp = strtok(NULL," ,'\"");
    if (tp != NULL){
      state[d] = atoi(tp);
      d++;
    }
  }
  if(d != DIM_STATE){
    fprintf(stderr, "The dimension of initial states is not adequate!\n");
    exit(1);
  }
}

int is_pruned(int *current, int *target, int length){
  // Calculate the number of steps required in each dimension and take the sum of them.
  double dist;
  for(int d = 0; d < DIM_STATE; d++){
    double comp = (double)(current[d]-target[d]);
    if(comp >= 0.0){
      double comp_pos = PRUNING_COEFF_POS[d] * comp;
      if(dist < comp_pos){
	dist = comp_pos;
      }
    }else{
      double comp_neg = PRUNING_COEFF_NEG[d] * comp;
      if(dist < comp_neg){
	dist = comp_neg;
      }
    }
  }
  // Judge whether the distance is over the length or not.
  if(dist <= length){
    return FALSE; // not to be pruned
  }else{
    return TRUE; // to be pruned
  }
}

//-----------------------------------------------------------------
// Euler
//-----------------------------------------------------------------
double evaluate_with_Euler(double dt, int max_M){
  int M;
  int64 num_of_terms, pre_num_of_terms;
  double working_coeffs;
  int *working_state;
  working_state = (int*)clear_alloc(sizeof(int), DIM_STATE);

  MAXIMUM_ID_NUMBER = 0;

  // Evalate constant rates (coefficients).
  for(int e = 0; e < NUM_ORI_EVENT; e++){
    double comp = ORI_RATE[e];
    for(int i = 0; i < NUM_CONSTANT; i++){
      comp = comp * pow(CONSTS[i], (double)ORI_INDEX4CONST[e][i]);
    }
    ORI_EVENT_CONST_RATE[e] = comp;
    printf("%f\n", ORI_EVENT_CONST_RATE[e]);
  }

  // Start from this M
  int initial_M = 2;
  M = initial_M;

  double pre_result = 0.0;
  double current_result;
  // Prepare the temporal output file for the sequences of differentials.
  FILE *fp;
  fp = f_open("log.txt", "w");
  fprintf(fp, "# T = %f\n", dt);
  fprintf(fp, "# [ ");
  for(int d = 0; d < DIM_STATE; d++){
    fprintf(fp, "%d ", INITIAL_STATE[d]);
  }
  fprintf(fp, "] -> [ ");
  for(int d = 0; d < DIM_STATE; d++){
    fprintf(fp, "%d ", TARGET_STATE[d]);
  }
  fprintf(fp, "]\n");
  fflush(fp);
  
  M=max_M; //######
  // while(M <= max_M){
    double dt_M = dt / (double)M;

    // Perform the initialization for each M.
    for(int d = 0; d < DIM_STATE; d++){
      STATE[0][d] = INITIAL_STATE[d];
    }
    COEFF[0] = 1.0;
    num_of_terms = 1;
    
    // Calculate coefficients and states iteratively.
    for(int depth = 1; depth <= M; depth++){
      // Copy the obtained states in the previous calculation.
      pre_num_of_terms = num_of_terms;
      for(int64 i = 0; i < pre_num_of_terms; i++){
	ORI_COEFF[i] = COEFF[i];
	for(int d = 0; d < DIM_STATE; d++){
	  ORI_STATE[i][d] = STATE[i][d];
	}
      }
      num_of_terms = 0;
      for(int64 i = 0; i < pre_num_of_terms; i++){
	for(int e1 = 0; e1 < NUM_EVENT; e1++){
	  for(int d = 0; d < DIM_STATE; d++){
	    working_state[d] = ORI_STATE[i][d];
	  }
	  working_coeffs = occur_event_Euler(e1, working_state, dt_M);	  
	  if(working_coeffs == 0.0){
	    continue;
	  }
	  
	  // Check whether the new state is a candidate or not.
	  if(is_pruned(working_state, TARGET_STATE, M-depth+1) == TRUE && IS_CUTOFF_USED == TRUE){
	    continue; // -------------------------------------------------- simple cutoff
	  }
	  
	  // Check whether the state has already appeared or not.
	  int is_already_appeared = FALSE;
	  int64 detected_index;
	  for(int64 j = 0; j < num_of_terms; j++){
	    int is_different = FALSE;
	    for(int d = 0; d < DIM_STATE; d++){
	      if(working_state[d] != STATE[j][d]){
		is_different = TRUE;
		break;
	      }
	    }
	    if(is_different == TRUE){
	      continue;
	    }else{
	      is_already_appeared = TRUE;
	      detected_index = j;
	      break;
	    }
	  }
	  if(is_already_appeared == TRUE){
	    // Revise the obtained coefficients.
	    // [CAUTION] This code could cause the problem of "cancellation of significant digits"...
	    COEFF[detected_index] = COEFF[detected_index] + ORI_COEFF[i]*working_coeffs;
      }else{
	    // Add the new index.
	    COEFF[num_of_terms] = ORI_COEFF[i]*working_coeffs;
	    for(int d = 0; d < DIM_STATE; d++){
	      STATE[num_of_terms][d] = working_state[d];
	    }
	    num_of_terms += 1;
	    if(num_of_terms >= MAX_ID){
	      fprintf(stderr, "ERROR: The number of IDs is insufficient. Use larger MAX_ID.\n");
	      exit(1);
	    }
	    if(num_of_terms >= MAXIMUM_ID_NUMBER){
	      MAXIMUM_ID_NUMBER = num_of_terms;
	    }
	  }
	} // end of event
      }
    }

    // printf("%f,%f,%f, %lld, %d, %d\n", COEFF[detected_index], ORI_COEFF[i], working_coeffs, i, STATE[i][0], STATE[i][1]);//#########

    // Check 
    for(int i = 0; i < num_of_terms; i++){
      // Check whether the state is the target one or not.
      int is_target_state = TRUE;
      for(int d = 0; d < DIM_STATE; d++){
	if(STATE[i][d] != TARGET_STATE[d]){
	  is_target_state = FALSE;
	  break;
	}
      }
      if(is_target_state == TRUE){
	current_result = COEFF[i];
	break;
      }
    }

    // Estimation
    double C, estimated_value;
    C = (double)(M*(M-1))/dt*(pre_result-current_result);
    estimated_value = current_result - C*dt/(double)M;
    
    fprintf(fp, "%3d %+.15e %.15e\n", M, current_result, estimated_value);
    fflush(fp);

    pre_result = current_result;
    M += 1;
  // } // end of one iteration for M
  M = M-1;
  
  fclose(fp);
  free(working_state);
  return current_result;
}


//-----------------------------------------------------------------
// Heun
//-----------------------------------------------------------------

double evaluate_with_Heun(double dt, int max_m){
  int M;
  int64 num_of_terms, pre_num_of_terms;
  double working_coeffs, working_coeffs_2nd;
  int *working_state, *working_state_2nd;
  working_state = (int*)clear_alloc(sizeof(int), DIM_STATE);
  working_state_2nd = (int*)clear_alloc(sizeof(int), DIM_STATE);

  MAXIMUM_ID_NUMBER = 0;

  // Evalate constant rates (coefficients).
  for(int e = 0; e < NUM_ORI_EVENT; e++){
    double comp = ORI_RATE[e];
    for(int i = 0; i < NUM_CONSTANT; i++){
      comp = comp * pow(CONSTS[i], (double)ORI_INDEX4CONST[e][i]);
    }
    ORI_EVENT_CONST_RATE[e] = comp;
  }

  // Start from this M
  int initial_M = 2;
  M = initial_M;

  double pre_result = 0.0;
  double current_result;
  // Prepare the temporal output file for the sequences of differentials.
  FILE *fp;
  fp = f_open("log.txt", "w");
  fprintf(fp, "# T = %f\n", dt);
  fprintf(fp, "# [ ");
  for(int d = 0; d < DIM_STATE; d++){
    fprintf(fp, "%d ", INITIAL_STATE[d]);
  }
  fprintf(fp, "] -> [ ");
  for(int d = 0; d < DIM_STATE; d++){
    fprintf(fp, "%d ", TARGET_STATE[d]);
  }
  fprintf(fp, "]\n");
  fflush(fp);
  
  while(M <= max_m){
    double dt_M = dt / (double)M;
    
    // Perform the initialization for each M.
    for(int d = 0; d < DIM_STATE; d++){
      STATE[0][d] = INITIAL_STATE[d];
    }
    COEFF[0] = 1.0;
    num_of_terms = 1;
    
    // Calculate coefficients and states iteratively.
    for(int depth = 1; depth <= M; depth++){

      // Copy the obtained states in the previous calculation.
      pre_num_of_terms = num_of_terms;
      for(int64 i = 0; i < pre_num_of_terms; i++){
	ORI_COEFF[i] = COEFF[i];
	for(int d = 0; d < DIM_STATE; d++){
	  ORI_STATE[i][d] = STATE[i][d];
	}
      }
      num_of_terms = 0;

      // Occur 1st event.
      for(int64 i = 0; i < pre_num_of_terms; i++){
	for(int e1 = 0; e1 < NUM_EVENT; e1++){
	  for(int d = 0; d < DIM_STATE; d++){
	    working_state[d] = ORI_STATE[i][d];
	  }
	  working_coeffs = occur_event_Euler(e1, working_state, dt_M);
	  
	  if(working_coeffs == 0.0){
	    continue;
	  }

	  // Check whether the new state is a candidate or not.
	  if(is_pruned(working_state, TARGET_STATE, 2*(M-depth+1)) == TRUE && IS_CUTOFF_USED == TRUE){
	    continue; // -------------------------------------------------- simple cutoff
	  }

	  // Check whether the state has already appeared or not.
	  int is_already_appeared = FALSE;
	  int64 detected_index;
	  for(int64 j = 0; j < num_of_terms; j++){
	    int is_different = FALSE;
	    for(int d = 0; d < DIM_STATE; d++){
	      if(working_state[d] != STATE[j][d]){
		is_different = TRUE;
		break;
	      }
	    }
	    if(is_different == TRUE){
	      continue;
	    }else{
	      is_already_appeared = TRUE;
	      detected_index = j;
	      break;
	    }
	  }
	  if(is_already_appeared == TRUE){
	    // Revise the obtained coefficients.
	    // [CAUTION] This code could cause the problem of "cancellation of significant digits"...
	    COEFF[detected_index] = COEFF[detected_index] + ORI_COEFF[i]*working_coeffs;
	  }else{
	    // Add the new index.
	    COEFF[num_of_terms] = ORI_COEFF[i]*working_coeffs;
	    for(int d = 0; d < DIM_STATE; d++){
	      STATE[num_of_terms][d] = working_state[d];
	    }
	    num_of_terms += 1;
	    if(num_of_terms >= MAX_ID){
	      fprintf(stderr, "ERROR: The number of IDs is insufficient. Use larger MAX_ID.\n");
	      exit(1);
	    }
	  }
	} // end of event 1

	// 2nd order process
	for(int e1 = 0; e1 < NUM_EVENT; e1++){
	  for(int d = 0; d < DIM_STATE; d++){
	    working_state[d] = ORI_STATE[i][d];
	  }
	  working_coeffs = occur_event_Heun_2nd_process(e1, working_state, dt_M);
	  if(working_coeffs == 0.0){
	    continue;
	  }
	  for(int e2 = 0; e2 < NUM_EVENT; e2++){
	    for(int d = 0; d < DIM_STATE; d++){
	      working_state_2nd[d] = working_state[d];
	    }
	    working_coeffs_2nd = occur_event_Heun_2nd_process(e2, working_state_2nd, 0.5*dt_M);
	    if(working_coeffs_2nd == 0.0){
	      continue;
	    }

	    // Check whether the new state is a candidate or not.
	    if(is_pruned(working_state, TARGET_STATE, 2*(M-depth+1)) == TRUE && IS_CUTOFF_USED == TRUE){
	      continue; // -------------------------------------------------- simple cutoff
	    }

	    // Check whether the state has already appeared or not.
	    int is_already_appeared = FALSE;
	    int64 detected_index;
	    for(int64 j = 0; j < num_of_terms; j++){
	      int is_different = FALSE;
	      for(int d = 0; d < DIM_STATE; d++){
		if(working_state_2nd[d] != STATE[j][d]){
		  is_different = TRUE;
		  break;
		}
	      }
	      if(is_different == TRUE){
		continue;
	      }else{
		is_already_appeared = TRUE;
		detected_index = j;
		break;
	      }
	    }
	    if(is_already_appeared == TRUE){
	      // Revise the obtained coefficients.
	      // [CAUTION] This code could cause the problem of "cancellation of significant digits"...
	      COEFF[detected_index] = COEFF[detected_index] + ORI_COEFF[i]*working_coeffs*working_coeffs_2nd;
	    }else{
	      // Add the new index.
	      COEFF[num_of_terms] = ORI_COEFF[i]*working_coeffs*working_coeffs_2nd;
	      for(int d = 0; d < DIM_STATE; d++){
		STATE[num_of_terms][d] = working_state_2nd[d];
	      }
	      num_of_terms += 1;
	      if(num_of_terms >= MAX_ID){
		fprintf(stderr, "ERROR: The number of IDs is insufficient. Use larger MAX_ID.\n");
		exit(1);
	      }
	      if(num_of_terms >= MAXIMUM_ID_NUMBER){
		MAXIMUM_ID_NUMBER = num_of_terms;
	      }
	    }
	  } // end of event 2
	} // end of event 1 for the 2nd order process

      }
      
    }
    
    // Check 
    for(int i = 0; i < num_of_terms; i++){
      // Check whether the state is the target one or not.
      int is_target_state = TRUE;
      for(int d = 0; d < DIM_STATE; d++){
	if(STATE[i][d] != TARGET_STATE[d]){
	  is_target_state = FALSE;
	  break;
	}
      }
      if(is_target_state == TRUE){
	current_result = COEFF[i];
	break;
      }
    }

    // Estimation
    double C, estimated_value;
    C = (double)(M*M*(M-1)*(M-1)) / ( dt*dt*((double)(M*M-(M-1)*(M-1))) ) * (pre_result-current_result);
    estimated_value = current_result - C*(dt/(double)M)*(dt/(double)M);
    
    fprintf(fp, "%3d %+.15e %.15e\n", M, current_result, estimated_value);
    fflush(fp);

    pre_result = current_result;
    M += 1;
  } // end of one iteration for M
  M = M-1;
  
  fclose(fp);
  free(working_state);
  return current_result;
}



double occur_event_Euler(int event, int *state, double dt_M){
  if(is_valid_next_state(event, state) == FALSE){
    return 0.0;
  }
  double rate = calculate_transition_matrix_element(event, state);
  if(event == ID_DIAGONAL_EVENT){
    return 1.0 + dt_M*rate;
  }else{
    for(int d = 0; d < DIM_STATE; d++){
      state[d] += EVENT[event][d];
    }
    return dt_M*rate;
  }
}

double occur_event_Heun_2nd_process(int event, int *state, double dt_M){
  if(is_valid_next_state(event, state) == FALSE){
    return 0.0;
  }
  double rate = calculate_transition_matrix_element(event, state);
  if(event == ID_DIAGONAL_EVENT){
    return dt_M*rate;
  }else{
    for(int d = 0; d < DIM_STATE; d++){
      state[d] += EVENT[event][d];
    }
    return dt_M*rate;
  }
}



//-----------------------------------------------------------------
// resolvent
//-----------------------------------------------------------------
double evaluate_with_resolvent(double dt, int max_m){
  int M;
  int64 num_of_terms, pre_num_of_terms;
  double working_coeffs;
  int *working_state;
  working_state = (int*)clear_alloc(sizeof(int), DIM_STATE);

  MAXIMUM_ID_NUMBER = 0;

  // Evalate constant rates (coefficients).
  for(int e = 0; e < NUM_ORI_EVENT; e++){
    double comp = ORI_RATE[e];
    for(int i = 0; i < NUM_CONSTANT; i++){
      comp = comp * pow(CONSTS[i], (double)ORI_INDEX4CONST[e][i]);
    }
    ORI_EVENT_CONST_RATE[e] = comp;
  }

  // Start from this M
  int initial_M = 2;
  M = initial_M;

  double pre_result = 0.0;
  double current_result;
  // Prepare the temporal output file for the sequences of differentials.
  FILE *fp;
  fp = f_open("log.txt", "w");
  fprintf(fp, "# T = %f\n", dt);
  fprintf(fp, "# [ ");
  for(int d = 0; d < DIM_STATE; d++){
    fprintf(fp, "%d ", INITIAL_STATE[d]);
  }
  fprintf(fp, "] -> [ ");
  for(int d = 0; d < DIM_STATE; d++){
    fprintf(fp, "%d ", TARGET_STATE[d]);
  }
  fprintf(fp, "]\n");
  fflush(fp);
  // M = 20;
  while(M <= max_m)
  {
    double dt_M = dt / (double)M;
    double first_rate, final_rate;
    
    // Calculate the first and final rates. These are multiplied in the final result.
    double rate = 0.0;
    for(int j = 0; j < LIST_LEN_E2ORI[ID_DIAGONAL_EVENT]; j++){
      int pre_e = LIST_START_INDEX_E2ORI[ID_DIAGONAL_EVENT] + j;
      double state_dependent_rate = 1.0;
      for(int d = 0; d < DIM_STATE; d++){
	for(int fact = 0; fact < ORI_INDEX4RATE[pre_e][d]; fact++){
	  state_dependent_rate = state_dependent_rate * (double)(INITIAL_STATE[d] - fact);
	}
      }
      rate = rate + state_dependent_rate * ORI_EVENT_CONST_RATE[pre_e];
    }
    first_rate = 1.0 - dt_M*rate;
    // Thr followings are not needed when the diagonal part does not include a constant.
    rate = 0.0;
    for(int j = 0; j < LIST_LEN_E2ORI[ID_DIAGONAL_EVENT]; j++){
      int pre_e = LIST_START_INDEX_E2ORI[ID_DIAGONAL_EVENT] + j;
      double state_dependent_rate = 1.0;
      for(int d = 0; d < DIM_STATE; d++){
	for(int fact = 0; fact < ORI_INDEX4RATE[pre_e][d]; fact++){
	  state_dependent_rate = state_dependent_rate * (double)(TARGET_STATE[d] - fact);
	}
      }
      rate = rate + state_dependent_rate * ORI_EVENT_CONST_RATE[pre_e];
    }
    final_rate = 1.0/(1.0 - dt_M*rate);
    
    // Perform the initialization for each M.
    for(int d = 0; d < DIM_STATE; d++){
      STATE[0][d] = INITIAL_STATE[d];
    }
    COEFF[0] = 1.0;
    num_of_terms = 1;
    
    // Calculate coefficients and states iteratively.
    for(int depth = 1; depth <= M; depth++){
      // Copy the obtained states in the previous calculation.
      pre_num_of_terms = num_of_terms;
      for(int64 i = 0; i < pre_num_of_terms; i++){
	ORI_COEFF[i] = COEFF[i];
	for(int d = 0; d < DIM_STATE; d++){
	  ORI_STATE[i][d] = STATE[i][d];
	}
      }
      num_of_terms = 0;
      for(int64 i = 0; i < pre_num_of_terms; i++){
	for(int e = 0; e < NUM_EVENT; e++){
	  for(int d = 0; d < DIM_STATE; d++){
	    working_state[d] = ORI_STATE[i][d];
	  }
	  working_coeffs = occur_event_resolvent(e, working_state, dt_M);
	  
	  if(working_coeffs == 0.0){
	    continue;
	  }
	  
	  // Check whether the new state is a candidate or not.
	  if(is_pruned(working_state, TARGET_STATE, M-depth+1) == TRUE && IS_CUTOFF_USED == TRUE){
	    continue; // -------------------------------------------------- simple cutoff
	  }
	  
	  // Check whether the state has already appeared or not.
	  int is_already_appeared = FALSE;
	  int64 detected_index;
	  for(int64 j = 0; j < num_of_terms; j++){
	    int is_different = FALSE;
	    for(int d = 0; d < DIM_STATE; d++){
	      if(working_state[d] != STATE[j][d]){
		is_different = TRUE;
		break;
	      }
	    }
	    if(is_different == TRUE){
	      continue;
	    }else{
	      is_already_appeared = TRUE;
	      detected_index = j;
	      break;
	    }
	  }
	  if(is_already_appeared == TRUE){
	    // Revise the obtained coefficients.
	    // [CAUTION] This code could cause the problem of "cancellation of significant digits"...
	    COEFF[detected_index] = COEFF[detected_index] + ORI_COEFF[i]*working_coeffs;
	  }else{
	    // Add the new index.
	    COEFF[num_of_terms] = ORI_COEFF[i]*working_coeffs;
	    for(int d = 0; d < DIM_STATE; d++){
	      STATE[num_of_terms][d] = working_state[d];
	    }
	    num_of_terms += 1;
	    if(num_of_terms >= MAX_ID){
	      fprintf(stderr, "ERROR: The number of IDs is insufficient. Use larger MAX_ID.\n");
	      exit(1);
	    }
	    if(num_of_terms >= MAXIMUM_ID_NUMBER){
	      MAXIMUM_ID_NUMBER = num_of_terms;
	    }
	  }
	}
      }
    }
    
    // Check 
    for(int i = 0; i < num_of_terms; i++){
      // Check whether the state is the target one or not.
      int is_target_state = TRUE;
      for(int d = 0; d < DIM_STATE; d++){
	if(STATE[i][d] != TARGET_STATE[d]){
	  is_target_state = FALSE;
	  break;
	}
      }
      if(is_target_state == TRUE){
	current_result = first_rate * COEFF[i] * final_rate;
	break;
      }
    }
    // Estimation
    double C, estimated_value;
    C = (double)(M*(M-1))/dt*(pre_result-current_result);
    estimated_value = current_result - C*dt/(double)M;
    
    fprintf(fp, "%3d %+.15e %.15e\n", M, current_result, estimated_value);
    fflush(fp);

    pre_result = current_result;
    M += 1;
  } // end of one iteration for M
  M = M-1;
  
  fclose(fp);
  free(working_state);
  return current_result;
}


double calculate_transition_matrix_element(int event, int *state){
  double rate = 0.0;
  for(int j = 0; j < LIST_LEN_E2ORI[event]; j++){
    int pre_e = LIST_START_INDEX_E2ORI[event] + j;
    double state_dependent_rate = 1.0;
    for(int d = 0; d < DIM_STATE; d++){
      for(int fact = 0; fact < ORI_INDEX4RATE[pre_e][d]; fact++){
	state_dependent_rate = state_dependent_rate * (double)(state[d] - fact);
      }
    }
    rate = rate + state_dependent_rate * ORI_EVENT_CONST_RATE[pre_e];
  }
  return rate;
}

int is_valid_next_state(int event, int *state){
  int is_valid_next_state = TRUE;
  for(int d = 0; d < DIM_STATE; d++){
    // If the next state is in the negative region, it is not valid.    
    if(state[d] + EVENT[event][d] < 0){
      is_valid_next_state = FALSE;
      break;
    }
  }
  return is_valid_next_state;
}


double occur_event_resolvent(int event, int *state, double dt_M){
  if(is_valid_next_state(event, state) == FALSE){
    return 0.0;
  }
  // Calculate the denominator.
  double rate;
  rate = calculate_transition_matrix_element(ID_DIAGONAL_EVENT, state);
  double denom_inv;
  denom_inv = 1.0/(1.0 - dt_M*rate);
  if(event == ID_DIAGONAL_EVENT){
    return denom_inv;
  }else{
    denom_inv = denom_inv*denom_inv;
    rate = calculate_transition_matrix_element(event, state);
    double numer = dt_M*rate;
    for(int d = 0; d < DIM_STATE; d++){
      state[d] += EVENT[event][d];
    }
    return numer * denom_inv;
  }  
  return rate;
}




//----------------------------------------------------------------------------------
//--- Crank Nicolson
//----------------------------------------------------------------------------------

void make_2nd_order_event(void){
  LIST_EVENTS_IN_DIAGONAL_RELAY_PATH = (int**)clear_alloc_2d(sizeof(int), NUM_EVENT, 2);
  LENGTH_OF_DIAGONAL_RELAY_PATH = 0;
  int len_list;
  int is_relay_path;
  for(int event = 0; event < NUM_EVENT; event++){
    if(event != ID_DIAGONAL_EVENT){ continue; }
    len_list = 0;
    for(int e1 = 0; e1 < NUM_EVENT; e1++){
      if(e1 == ID_DIAGONAL_EVENT){ continue; }
      for(int e2 = 0; e2 < NUM_EVENT; e2++){
	if(e2 == ID_DIAGONAL_EVENT){ continue; }
	is_relay_path = TRUE;
	for(int d = 0; d < DIM_STATE; d++){
	  if(EVENT[event][d] != (EVENT[e1][d] + EVENT[e2][d])){
	    is_relay_path = FALSE;
	    break;
	  }
	}
	if(is_relay_path == FALSE){ continue; }
	// e1 -> e2 is a relay path
	LIST_EVENTS_IN_DIAGONAL_RELAY_PATH[len_list][0] = e1;
	LIST_EVENTS_IN_DIAGONAL_RELAY_PATH[len_list][1] = e2;
	len_list++;
      }
    }
    LENGTH_OF_DIAGONAL_RELAY_PATH += len_list;
  }
}

void free_memories_allocated_in_make_2nd_order_event(void){
  free_2d(LIST_EVENTS_IN_DIAGONAL_RELAY_PATH);
}


double evaluate_with_Crank_Nicolson_simple(double dt, int max_m){
  int M;
  int64 num_of_terms, pre_num_of_terms;
  double working_coeffs, working_coeffs_2nd;
  int *working_state, *working_state_2nd;
  working_state = (int*)clear_alloc(sizeof(int), DIM_STATE);
  working_state_2nd = (int*)clear_alloc(sizeof(int), DIM_STATE);

  MAXIMUM_ID_NUMBER = 0;

  // Evalate constant rates (coefficients).
  for(int e = 0; e < NUM_ORI_EVENT; e++){
    double comp = ORI_RATE[e];
    for(int i = 0; i < NUM_CONSTANT; i++){
      comp = comp * pow(CONSTS[i], (double)ORI_INDEX4CONST[e][i]);
    }
    ORI_EVENT_CONST_RATE[e] = comp;
  }

  // Start from this M
  int initial_M = 2;
  M = initial_M;

  double pre_result = 0.0;
  double current_result;
  // Prepare the temporal output file for the sequences of differentials.
  FILE *fp;
  fp = f_open("log.txt", "w");
  fprintf(fp, "# T = %f\n", dt);
  fprintf(fp, "# [ ");
  for(int d = 0; d < DIM_STATE; d++){
    fprintf(fp, "%d ", INITIAL_STATE[d]);
  }
  fprintf(fp, "] -> [ ");
  for(int d = 0; d < DIM_STATE; d++){
    fprintf(fp, "%d ", TARGET_STATE[d]);
  }
  fprintf(fp, "]\n");
  fflush(fp);

  
  while(M <= max_m){
    double dt_M = dt / (double)M;
    
    // Perform the initialization for each M.
    for(int d = 0; d < DIM_STATE; d++){
      STATE[0][d] = INITIAL_STATE[d];
    }
    COEFF[0] = 1.0;
    num_of_terms = 1;
    
    // Calculate coefficients and states iteratively.
    for(int depth = 1; depth <= M; depth++){

      // Copy the obtained states in the previous calculation.
      pre_num_of_terms = num_of_terms;
      for(int64 i = 0; i < pre_num_of_terms; i++){
	ORI_COEFF[i] = COEFF[i];
	for(int d = 0; d < DIM_STATE; d++){
	  ORI_STATE[i][d] = STATE[i][d];
	}
      }
      num_of_terms = 0;

      { // start of 1st step for the Clank-Nicolson
	for(int64 i = 0; i < pre_num_of_terms; i++){
	  for(int e = 0; e < NUM_EVENT; e++){
	    for(int d = 0; d < DIM_STATE; d++){
	      working_state[d] = ORI_STATE[i][d];
	    }
	    working_coeffs = occur_event_Euler(e, working_state, 0.5*dt_M);
	    if(working_coeffs == 0.0){
	      continue;
	    }
	    
	    // Check whether the state has already appeared or not.
	    int is_already_appeared = FALSE;
	    int64 detected_index;
	    for(int64 j = 0; j < num_of_terms; j++){
	      int is_different = FALSE;
	      for(int d = 0; d < DIM_STATE; d++){
		if(working_state[d] != STATE[j][d]){
		  is_different = TRUE;
		  break;
		}
	      }
	      if(is_different == TRUE){
		continue;
	      }else{
		is_already_appeared = TRUE;
		detected_index = j;
		break;
	      }
	    }
	    if(is_already_appeared == TRUE){
	      // Revise the obtained coefficients.
	      // [CAUTION] This code could cause the problem of "cancellation of significant digits"...
	      COEFF[detected_index] = COEFF[detected_index] + ORI_COEFF[i]*working_coeffs;
	    }else{
	      // Add the new index.
	      COEFF[num_of_terms] = ORI_COEFF[i]*working_coeffs;
	      for(int d = 0; d < DIM_STATE; d++){
		STATE[num_of_terms][d] = working_state[d];
	      }
	      num_of_terms += 1;
	      if(num_of_terms >= MAX_ID){
		fprintf(stderr, "ERROR: The number of IDs is insufficient. Use larger MAX_ID.\n");
		exit(1);
	      }
	    }
	  }
	}
      } // end of 1st step for the Crank-Nicolson

      
      // Copy the obtained states in the previous calculation.
      pre_num_of_terms = num_of_terms;
      for(int64 i = 0; i < pre_num_of_terms; i++){
	ORI_COEFF[i] = COEFF[i];
	for(int d = 0; d < DIM_STATE; d++){
	  ORI_STATE[i][d] = STATE[i][d];
	}
      }
      num_of_terms = 0;
      
      { // start of 2nd step for the Clank-Nicolson
	for(int64 i = 0; i < pre_num_of_terms; i++){
	  for(int e1 = 0; e1 < NUM_EVENT; e1++){
	    for(int d = 0; d < DIM_STATE; d++){
	      working_state[d] = ORI_STATE[i][d];
	    }
	    working_coeffs = occur_event_resolvent_2nd_order(e1, working_state, 0.5*dt_M);
	    
	    if(working_coeffs == 0.0){
	      continue;
	    }

	    // Check whether the new state is a candidate or not.
	    if(is_pruned(working_state, TARGET_STATE, 3*(M-depth+1)) == TRUE && IS_CUTOFF_USED == TRUE){
	      continue; // -------------------------------------------------- simple cutoff
	    }

	    // Check whether the state has already appeared or not.
	    int is_already_appeared = FALSE;
	    int64 detected_index;
	    for(int64 j = 0; j < num_of_terms; j++){
	      int is_different = FALSE;
	      for(int d = 0; d < DIM_STATE; d++){
		if(working_state[d] != STATE[j][d]){
		  is_different = TRUE;
		  break;
		}
	      }
	      if(is_different == TRUE){
		continue;
	      }else{
		is_already_appeared = TRUE;
		detected_index = j;
		break;
	      }
	    }
	    if(is_already_appeared == TRUE){
	      // Revise the obtained coefficients.
	      // [CAUTION] This code could cause the problem of "cancellation of significant digits"...
	      COEFF[detected_index] = COEFF[detected_index] + ORI_COEFF[i]*working_coeffs;
	    }else{
	      // Add the new index.
	      COEFF[num_of_terms] = ORI_COEFF[i]*working_coeffs;
	      for(int d = 0; d < DIM_STATE; d++){
		STATE[num_of_terms][d] = working_state[d];
	      }
	      num_of_terms += 1;
	      if(num_of_terms >= MAX_ID){
		fprintf(stderr, "ERROR: The number of IDs is insufficient. Use larger MAX_ID.\n");
		exit(1);
	      }
	    }
	  } // end of event 1
	  
	  // 2nd order process
	  for(int e1 = 0; e1 < NUM_EVENT; e1++){
	    if(e1 == ID_DIAGONAL_EVENT){ continue; }
	    for(int d = 0; d < DIM_STATE; d++){
	      working_state[d] = ORI_STATE[i][d];
	    }
	    working_coeffs = occur_event_Euler(e1, working_state, 0.5*dt_M);
	    if(working_coeffs == 0.0){
	      continue;
	    }
	    for(int e2 = 0; e2 < NUM_EVENT; e2++){
	      if(e2 == ID_DIAGONAL_EVENT){ continue; }
	      for(int d = 0; d < DIM_STATE; d++){
		working_state_2nd[d] = working_state[d];
	      }
	      
	      int is_diagonal = TRUE;
	      for(int d = 0; d < DIM_STATE; d++){
		if(working_state_2nd[d]+EVENT[e2][d] != ORI_STATE[i][d]){
		  is_diagonal = FALSE;
		  break;
		}
	      }
	      if(is_diagonal == TRUE){
		continue;
	      }
	      	      
	      working_coeffs_2nd = occur_event_Euler(e2, working_state_2nd, 0.5*dt_M);
	      if(working_coeffs_2nd == 0.0){
		continue;
	      }

	      // Check whether the new state is a candidate or not.
	      if(is_pruned(working_state, TARGET_STATE, 3*(M-depth+1)) == TRUE && IS_CUTOFF_USED == TRUE){
		continue; // -------------------------------------------------- simple cutoff
	      }

	      // Check whether the state has already appeared or not.
	      int is_already_appeared = FALSE;
	      int64 detected_index;
	      for(int64 j = 0; j < num_of_terms; j++){
		int is_different = FALSE;
		for(int d = 0; d < DIM_STATE; d++){
		  if(working_state_2nd[d] != STATE[j][d]){
		    is_different = TRUE;
		    break;
		  }
		}
		if(is_different == TRUE){
		  continue;
		}else{
		  is_already_appeared = TRUE;
		  detected_index = j;
		  break;
		}
	      }
	      if(is_already_appeared == TRUE){
		// Revise the obtained coefficients.
		// [CAUTION] This code could cause the problem of "cancellation of significant digits"...
		COEFF[detected_index] = COEFF[detected_index] + ORI_COEFF[i]*working_coeffs*working_coeffs_2nd;
	      }else{
		// Add the new index.
		COEFF[num_of_terms] = ORI_COEFF[i]*working_coeffs*working_coeffs_2nd;
		for(int d = 0; d < DIM_STATE; d++){
		  STATE[num_of_terms][d] = working_state_2nd[d];
		}
		num_of_terms += 1;
		if(num_of_terms >= MAX_ID){
		  fprintf(stderr, "ERROR: The number of IDs is insufficient. Use larger MAX_ID.\n");
		  exit(1);
		}
		if(num_of_terms >= MAXIMUM_ID_NUMBER){
		  MAXIMUM_ID_NUMBER = num_of_terms;
		}
		
	      }	      
	    } // end of event 2
	  } // end of event 1 for the 2nd order process
	}
      } // end of 2nd step
      
    }
    
    // Check 
    for(int i = 0; i < num_of_terms; i++){
      // Check whether the state is the target one or not.
      int is_target_state = TRUE;
      for(int d = 0; d < DIM_STATE; d++){
	if(STATE[i][d] != TARGET_STATE[d]){
	  is_target_state = FALSE;
	  break;
	}
      }
      if(is_target_state == TRUE){
	current_result = COEFF[i];
	break;
      }
    }
    
    // Estimation
    double C, estimated_value;
    C = (double)(M*M*(M-1)*(M-1)) / ( dt*dt*((double)(M*M-(M-1)*(M-1))) ) * (pre_result-current_result);
    estimated_value = current_result - C*(dt/(double)M)*(dt/(double)M);
    
    fprintf(fp, "%3d %+.15e %.15e\n", M, current_result, estimated_value);
    fflush(fp);

    pre_result = current_result;
    M += 1;
  } // end of one iteration for M
  M = M-1;
  
  fclose(fp);
  free(working_state);
  free(working_state_2nd);
  return current_result;
}








double occur_event_resolvent_2nd_order(int event, int *state, double dt_M){
  if(is_valid_next_state(event, state) == FALSE){
    return 0.0;
  }

  // Diagonal case
  if(event == ID_DIAGONAL_EVENT){
    double final_rate;
    double rate_1st_order = calculate_transition_matrix_element(ID_DIAGONAL_EVENT, state);
    double rate_2nd_order = 0.0;
    for(int l = 0; l < LENGTH_OF_DIAGONAL_RELAY_PATH; l++){
      double rate_1, rate_2;
      int relay_event_1 = LIST_EVENTS_IN_DIAGONAL_RELAY_PATH[l][0];
      int relay_event_2 = LIST_EVENTS_IN_DIAGONAL_RELAY_PATH[l][1];
      if(is_valid_next_state(relay_event_1, state) == FALSE){ continue; }
      rate_1 = calculate_transition_matrix_element(relay_event_1, state);
      for(int d = 0; d < DIM_STATE; d++){
	state[d] += EVENT[relay_event_1][d];
      }
      rate_2 = calculate_transition_matrix_element(relay_event_2, state);
      for(int d = 0; d < DIM_STATE; d++){
	state[d] -= EVENT[relay_event_1][d];
      }
      rate_2nd_order += rate_1*rate_2;
    }
    final_rate = 1.0/(1.0 - dt_M*rate_1st_order - dt_M*dt_M*rate_2nd_order);
    return final_rate;
  }

  // Non-diagonal case
  double final_rate;
  double denom_before;
  double denom_after;
  {
    // Evaluate the denominator before transition
    double rate_1st_order = calculate_transition_matrix_element(ID_DIAGONAL_EVENT, state);
    denom_before = 1.0/(1.0 - dt_M*rate_1st_order);
  }
  {
    // Evaluate the denominator after transition
    //-- move to the state after the event
    for(int d = 0; d < DIM_STATE; d++){
      state[d] += EVENT[event][d];
    }
    double rate_1st_order = calculate_transition_matrix_element(ID_DIAGONAL_EVENT, state);
    denom_after = 1.0/(1.0 - dt_M*rate_1st_order);
    //-- move back to the original state
    for(int d = 0; d < DIM_STATE; d++){
      state[d] -= EVENT[event][d];
    }
  }
  
  double rate_1st_order = calculate_transition_matrix_element(event, state);
  final_rate = dt_M*rate_1st_order*denom_before*denom_after;
  for(int d = 0; d < DIM_STATE; d++){
    state[d] += EVENT[event][d];
  }
  return final_rate;
}







//---------------------------
double evaluate_with_Crank_Nicolson(double dt, int max_m){
  int M;
  int64 num_of_terms, pre_num_of_terms;
  double working_coeffs, working_coeffs_2nd;
  int *working_state, *working_state_2nd;
  working_state = (int*)clear_alloc(sizeof(int), DIM_STATE);
  working_state_2nd = (int*)clear_alloc(sizeof(int), DIM_STATE);

  MAXIMUM_ID_NUMBER = 0;

  // Evalate constant rates (coefficients).
  for(int e = 0; e < NUM_ORI_EVENT; e++){
    double comp = ORI_RATE[e];
    for(int i = 0; i < NUM_CONSTANT; i++){
      comp = comp * pow(CONSTS[i], (double)ORI_INDEX4CONST[e][i]);
    }
    ORI_EVENT_CONST_RATE[e] = comp;
  }

  // Start from this M
  int initial_M = 2;
  M = initial_M;

  double pre_result = 0.0;
  double current_result;
  // Prepare the temporal output file for the sequences of differentials.
  FILE *fp;
  fp = f_open("log.txt", "w");
  fprintf(fp, "# T = %f\n", dt);
  fprintf(fp, "# [ ");
  for(int d = 0; d < DIM_STATE; d++){
    fprintf(fp, "%d ", INITIAL_STATE[d]);
  }
  fprintf(fp, "] -> [ ");
  for(int d = 0; d < DIM_STATE; d++){
    fprintf(fp, "%d ", TARGET_STATE[d]);
  }
  fprintf(fp, "]\n");
  fflush(fp);

  
  while(M <= max_m){
    double dt_M = dt / (double)M;
    
    // Perform the initialization for each M.
    for(int d = 0; d < DIM_STATE; d++){
      STATE[0][d] = INITIAL_STATE[d];
    }
    COEFF[0] = 1.0;
    num_of_terms = 1;
    
    // Calculate coefficients and states iteratively.
    for(int depth = 1; depth <= M; depth++){

      // Copy the obtained states in the previous calculation.
      pre_num_of_terms = num_of_terms;
      for(int64 i = 0; i < pre_num_of_terms; i++){
	ORI_COEFF[i] = COEFF[i];
	for(int d = 0; d < DIM_STATE; d++){
	  ORI_STATE[i][d] = STATE[i][d];
	}
      }
      num_of_terms = 0;

      { // start of 1st step for the Clank-Nicolson
	for(int64 i = 0; i < pre_num_of_terms; i++){
	  for(int e = 0; e < NUM_EVENT; e++){
	    for(int d = 0; d < DIM_STATE; d++){
	      working_state[d] = ORI_STATE[i][d];
	    }
	    working_coeffs = occur_event_Euler(e, working_state, 0.5*dt_M);
	    if(working_coeffs == 0.0){
	      continue;
	    }
	    
	    // Check whether the state has already appeared or not.
	    int is_already_appeared = FALSE;
	    int64 detected_index;
	    for(int64 j = 0; j < num_of_terms; j++){
	      int is_different = FALSE;
	      for(int d = 0; d < DIM_STATE; d++){
		if(working_state[d] != STATE[j][d]){
		  is_different = TRUE;
		  break;
		}
	      }
	      if(is_different == TRUE){
		continue;
	      }else{
		is_already_appeared = TRUE;
		detected_index = j;
		break;
	      }
	    }
	    if(is_already_appeared == TRUE){
	      // Revise the obtained coefficients.
	      // [CAUTION] This code could cause the problem of "cancellation of significant digits"...
	      COEFF[detected_index] = COEFF[detected_index] + ORI_COEFF[i]*working_coeffs;
	    }else{
	      // Add the new index.
	      COEFF[num_of_terms] = ORI_COEFF[i]*working_coeffs;
	      for(int d = 0; d < DIM_STATE; d++){
		STATE[num_of_terms][d] = working_state[d];
	      }
	      num_of_terms += 1;
	      if(num_of_terms >= MAX_ID){
		fprintf(stderr, "ERROR: The number of IDs is insufficient. Use larger MAX_ID.\n");
		exit(1);
	      }
	    }
	  }
	}
      } // end of 1st step for the Crank-Nicolson

      
      // Copy the obtained states in the previous calculation.
      pre_num_of_terms = num_of_terms;
      for(int64 i = 0; i < pre_num_of_terms; i++){
	ORI_COEFF[i] = COEFF[i];
	for(int d = 0; d < DIM_STATE; d++){
	  ORI_STATE[i][d] = STATE[i][d];
	}
      }
      num_of_terms = 0;
      
      { // start of 2nd step for the Clank-Nicolson
	for(int64 i = 0; i < pre_num_of_terms; i++){
	  for(int e1 = 0; e1 < NUM_EVENT; e1++){
	    for(int d = 0; d < DIM_STATE; d++){
	      working_state[d] = ORI_STATE[i][d];
	    }
	    working_coeffs = occur_event_resolvent_2nd_order(e1, working_state, 0.5*dt_M);
	    
	    if(working_coeffs == 0.0){
	      continue;
	    }

	    // Check whether the new state is a candidate or not.
	    if(is_pruned(working_state, TARGET_STATE, 3*(M-depth+1)) == TRUE && IS_CUTOFF_USED == TRUE){
	      continue; // -------------------------------------------------- simple cutoff
	    }

	    // Check whether the state has already appeared or not.
	    int is_already_appeared = FALSE;
	    int64 detected_index;
	    for(int64 j = 0; j < num_of_terms; j++){
	      int is_different = FALSE;
	      for(int d = 0; d < DIM_STATE; d++){
		if(working_state[d] != STATE[j][d]){
		  is_different = TRUE;
		  break;
		}
	      }
	      if(is_different == TRUE){
		continue;
	      }else{
		is_already_appeared = TRUE;
		detected_index = j;
		break;
	      }
	    }
	    if(is_already_appeared == TRUE){
	      // Revise the obtained coefficients.
	      // [CAUTION] This code could cause the problem of "cancellation of significant digits"...
	      COEFF[detected_index] = COEFF[detected_index] + ORI_COEFF[i]*working_coeffs;
	    }else{
	      // Add the new index.
	      COEFF[num_of_terms] = ORI_COEFF[i]*working_coeffs;
	      for(int d = 0; d < DIM_STATE; d++){
		STATE[num_of_terms][d] = working_state[d];
	      }
	      num_of_terms += 1;
	      if(num_of_terms >= MAX_ID){
		fprintf(stderr, "ERROR: The number of IDs is insufficient. Use larger MAX_ID.\n");
		exit(1);
	      }
	    }
	  } // end of event 1
	  
	  // 2nd order process
	  double factor1, factor2, factor3, factor;
	  for(int e1 = 0; e1 < NUM_EVENT; e1++){
	    if(e1 == ID_DIAGONAL_EVENT){ continue; }
	    for(int d = 0; d < DIM_STATE; d++){
	      working_state[d] = ORI_STATE[i][d];
	    }
	    factor1 = evaluate_diagonal_factor(working_state, 0.5*dt_M);
	    working_coeffs = occur_event_Euler(e1, working_state, 0.5*dt_M);
	    if(working_coeffs == 0.0){
	      continue;
	    }
	    for(int e2 = 0; e2 < NUM_EVENT; e2++){
	      if(e2 == ID_DIAGONAL_EVENT){ continue; }
	      for(int d = 0; d < DIM_STATE; d++){
		working_state_2nd[d] = working_state[d];
	      }
	      
	      int is_diagonal = TRUE;
	      for(int d = 0; d < DIM_STATE; d++){
		if(working_state_2nd[d]+EVENT[e2][d] != ORI_STATE[i][d]){
		  is_diagonal = FALSE;
		  break;
		}
	      }
	      if(is_diagonal == TRUE){
		continue;
	      }

	      factor2 = evaluate_diagonal_factor(working_state_2nd, 0.5*dt_M);
	      working_coeffs_2nd = occur_event_Euler(e2, working_state_2nd, 0.5*dt_M);
	      if(working_coeffs_2nd == 0.0){
		continue;
	      }
	      factor3 = evaluate_diagonal_factor(working_state_2nd, 0.5*dt_M);

	      // Check whether the new state is a candidate or not.
	      if(is_pruned(working_state, TARGET_STATE, 3*(M-depth+1)) == TRUE && IS_CUTOFF_USED == TRUE){
		continue; // -------------------------------------------------- simple cutoff
	      }

	      factor = factor1*factor2*factor3;
	      
	      // Check whether the state has already appeared or not.
	      int is_already_appeared = FALSE;
	      int64 detected_index;
	      for(int64 j = 0; j < num_of_terms; j++){
		int is_different = FALSE;
		for(int d = 0; d < DIM_STATE; d++){
		  if(working_state_2nd[d] != STATE[j][d]){
		    is_different = TRUE;
		    break;
		  }
		}
		if(is_different == TRUE){
		  continue;
		}else{
		  is_already_appeared = TRUE;
		  detected_index = j;
		  break;
		}
	      }
	      if(is_already_appeared == TRUE){
		// Revise the obtained coefficients.
		// [CAUTION] This code could cause the problem of "cancellation of significant digits"...
		COEFF[detected_index] = COEFF[detected_index] + ORI_COEFF[i]*working_coeffs*working_coeffs_2nd*factor;
	      }else{
		// Add the new index.
		COEFF[num_of_terms] = ORI_COEFF[i]*working_coeffs*working_coeffs_2nd*factor;
		for(int d = 0; d < DIM_STATE; d++){
		  STATE[num_of_terms][d] = working_state_2nd[d];
		}
		num_of_terms += 1;
		if(num_of_terms >= MAX_ID){
		  fprintf(stderr, "ERROR: The number of IDs is insufficient. Use larger MAX_ID.\n");
		  exit(1);
		}
		if(num_of_terms >= MAXIMUM_ID_NUMBER){
		  MAXIMUM_ID_NUMBER = num_of_terms;
		}
		
	      }	      
	    } // end of event 2
	  } // end of event 1 for the 2nd order process
	}
      } // end of 2nd step
      
    }
    
    // Check 
    for(int i = 0; i < num_of_terms; i++){
      // Check whether the state is the target one or not.
      int is_target_state = TRUE;
      for(int d = 0; d < DIM_STATE; d++){
	if(STATE[i][d] != TARGET_STATE[d]){
	  is_target_state = FALSE;
	  break;
	}
      }
      if(is_target_state == TRUE){
	current_result = COEFF[i];
	break;
      }
    }
    
    // Estimation
    double C, estimated_value;
    C = (double)(M*M*(M-1)*(M-1)) / ( dt*dt*((double)(M*M-(M-1)*(M-1))) ) * (pre_result-current_result);
    estimated_value = current_result - C*(dt/(double)M)*(dt/(double)M);
    
    fprintf(fp, "%3d %+.15e %.15e\n", M, current_result, estimated_value);
    fflush(fp);

    pre_result = current_result;
    M += 1;
  } // end of one iteration for M
  M = M-1;
  
  fclose(fp);
  free(working_state);
  free(working_state_2nd);
  return current_result;
}




double evaluate_diagonal_factor(int *state, double dt_M){
  return 1.0/(1.0 - dt_M*calculate_transition_matrix_element(ID_DIAGONAL_EVENT, state));
}



// -----------------------------------------------------------
// Subroutines for file open and memory allocation
// -----------------------------------------------------------
FILE *f_open(const char *path, const char *mode){
  static FILE *pf_ret;
  pf_ret = fopen(path, mode);
  if( pf_ret == NULL ){
    fprintf( stderr, "can't open a file, %s.\n", path );
    exit( EXIT_FAILURE );
  }
  return pf_ret;
}

void *clear_alloc(size_t size, size_t nobj){
  void *pret;
  pret = calloc(nobj, size);
  if (pret == NULL) {
    fprintf(stderr,"memory allocation failure!");
    exit(0);
  }
  return pret;
}

void **clear_alloc_2d(size_t size, size_t row, size_t column){
  void **ppret = NULL;
  size_t i;
  ppret = (void **) clear_alloc(sizeof(void *), row);
  ppret[0] = (void *) clear_alloc(size, row * column);
  for (i = 0; i < row; i++) {
    ppret[i] = (void *) ((size_t) ppret[0] + i * size * column);
  }
  return ppret;
}

