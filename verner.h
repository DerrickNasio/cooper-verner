/*
 * verner.h
 */
#define ORDER 8
#define NUMBER_OF_STAGES 11
#define sqrt21 4.582575694955840006588047193728

/* State struct definition */
typedef struct
{
    double *k[NUMBER_OF_STAGES];
    double *y_temp;
} verner_state_t;

/* Function prototypes */
void verner(
    size_t number_of_equations,
    double t,
    double h_step,
    double y[number_of_equations],
    void (*RHS)(double, double *, double *),
    size_t number_of_outputs,
    double t_output[number_of_outputs]);

static void *
verner_alloc(size_t number_of_equations);
static void
verner_compute_stages(void *vstate,
                      size_t number_of_equations,
                      double t,
                      double h_step,
                      double y[number_of_equations],
                      void (*RHS)(double, double *, double *));
static void
verner_dense_output(void *vstate,
                    size_t number_of_equations,
                    double t,
                    double theta,
                    double h_step,
                    double y[number_of_equations]);
static void
verner_apply(void *vstate,
             size_t number_of_equations,
             double h_step,
             double y[number_of_equations]);
static void
verner_free(void *vstate);

void timestamp(void);
