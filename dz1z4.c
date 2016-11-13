
// bmp.h

#include "stdio.h"
#include "stdlib.h"

typedef struct{
    unsigned char B;
    unsigned char G;
    unsigned char R;
} RGB;

typedef struct {
    unsigned int filesz;
    unsigned short creator1;
    unsigned short creator2;
    unsigned int bmp_offset;
} bmpfile_header_t;

typedef struct {
    unsigned int header_sz;
    unsigned int width;
    unsigned int height;
    unsigned short nplanes;
    unsigned short bitspp;
    unsigned int compress_type;
    unsigned int bmp_bytesz;
    unsigned int hres;
    unsigned int vres;
    unsigned int ncolors;
    unsigned int nimpcolors;
} bmp_dib_header_t;

typedef enum {
    BI_RGB = 0,
    BI_RLE8,
    BI_RLE4,
    BI_BITFIELDS,
    BI_JPEG,
    BI_PNG,
} bmp_compression_method_t;

typedef struct{
    unsigned char magic[2];
    bmpfile_header_t file_header;
    bmp_dib_header_t dib_header;
    unsigned int* palette;
    void* pixel_map;
} bmp_image;

void create_bmp(RGB* bitmap, int height, int width, const char* filename){
    bmp_image image;

    int padded_width = 4*(((width*24)+31)/32);
    padded_width -= width*sizeof(RGB);

    char* pad = (char*) calloc (padded_width, sizeof(char));

    image.magic[0]='B';
    image.magic[1]='M';

    image.file_header.filesz = 2*sizeof(char) + sizeof(bmpfile_header_t) + sizeof(bmp_dib_header_t) + height*width*sizeof(RGB);
    image.file_header.creator1 = image.file_header.creator2 = 0;
    image.file_header.bmp_offset = 2*sizeof(char) + sizeof(bmpfile_header_t) + sizeof(bmp_dib_header_t);

    image.dib_header.header_sz = 40;//sizeof(bmp_dib_header_t);
    image.dib_header.width = width;
    image.dib_header.height = height;
    image.dib_header.nplanes = 1;
    image.dib_header.bitspp = 24;
    image.dib_header.compress_type = 0;
    image.dib_header.bmp_bytesz = width*height*sizeof(RGB);
    image.dib_header.hres = 0;
    image.dib_header.vres = 0;
    image.dib_header.ncolors = 0;
    image.dib_header.nimpcolors = 0;

    FILE* out_file = fopen(filename,"wb");

    fwrite(image.magic,sizeof(char),2,out_file);
    fwrite(&(image.file_header),sizeof(char),sizeof(bmpfile_header_t),out_file);
    fwrite(&(image.dib_header),sizeof(char),sizeof(bmp_dib_header_t),out_file);

    int h;
    for (h = height-1; h >= 0; h--){
        fwrite(&bitmap[h*width],sizeof(RGB),width,out_file);
        fwrite(pad,sizeof(char),padded_width,out_file);
    }

    fclose(out_file);
}

// end bmp.h


// utils.h

#ifndef _HEADER
#define _HEADER

#ifdef __cplusplus
extern "C" {
#endif

#include <unistd.h>

/* Command line parameters for benchmarks */
struct pb_Parameters {
    char *outFile;		/* If not NULL, the raw output of the
				 * computation should be saved to this
				 * file. The string is owned. */
    char **inpFiles;		/* A NULL-terminated array of strings
				 * holding the input file(s) for the
				 * computation.  The array and strings
				 * are owned. */
};

/* Read command-line parameters.
 *
 * The argc and argv parameters to main are read, and any parameters
 * interpreted by this function are removed from the argument list.
 *
 * A new instance of struct pb_Parameters is returned.
 * If there is an error, then an error message is printed on stderr
 * and NULL is returned.
 */
struct pb_Parameters *
pb_ReadParameters(int *_argc, char **argv);

/* Free an instance of struct pb_Parameters.
 */
void
pb_FreeParameters(struct pb_Parameters *p);

/* Count the number of input files in a pb_Parameters instance.
 */
int
pb_Parameters_CountInputs(struct pb_Parameters *p);

/* A time or duration. */
#if _POSIX_VERSION >= 200112L
typedef unsigned long long pb_Timestamp; /* time in microseconds */
#else
# error "Timestamps not implemented"
#endif

enum pb_TimerState {
    pb_Timer_STOPPED,
    pb_Timer_RUNNING,
};

struct pb_Timer {
    enum pb_TimerState state;
    pb_Timestamp elapsed;		/* Amount of time elapsed so far */
    pb_Timestamp init;		/* Beginning of the current time interval,
				 * if state is RUNNING.  End of the last
				 * recorded time interfal otherwise.  */
};

/* Reset a timer.
 * Use this to initialize a timer or to clear
 * its elapsed time.  The reset timer is stopped.
 */
void
pb_ResetTimer(struct pb_Timer *timer);

/* Start a timer.  The timer is set to RUNNING mode and
 * time elapsed while the timer is running is added to
 * the timer.
 * The timer should not already be running.
 */
void
pb_StartTimer(struct pb_Timer *timer);

/* Stop a timer.
 * This stops adding elapsed time to the timer.
 * The timer should not already be stopped.
 */
void
pb_StopTimer(struct pb_Timer *timer);

/* Get the elapsed time in seconds. */
double
pb_GetElapsedTime(struct pb_Timer *timer);

/* Execution time is assigned to one of these categories. */
enum pb_TimerID {
    pb_TimerID_NONE = 0,
    pb_TimerID_IO,		/* Time spent in input/output */
    pb_TimerID_KERNEL,		/* Time spent computing on the device,
				 * recorded asynchronously */
    pb_TimerID_COPY,		/* Time spent synchronously moving data
				 * to/from device and allocating/freeing
				 * memory on the device */
    pb_TimerID_DRIVER,		/* Time spent in the host interacting with the
				 * driver, primarily for recording the time
                                 * spent queueing asynchronous operations */
    pb_TimerID_COPY_ASYNC,	/* Time spent in asynchronous transfers */
    pb_TimerID_COMPUTE,		/* Time for all program execution other
				 * than parsing command line arguments,
				 * I/O, kernel, and copy */
    pb_TimerID_OVERLAP,		/* Time double-counted in asynchronous and
				 * host activity: automatically filled in,
				 * not intended for direct usage */
    pb_TimerID_LAST		/* Number of timer IDs */
};

/* Dynamic list of asynchronously tracked times between events */
struct pb_async_time_marker_list {
    char *label; // actually just a pointer to a string
    enum pb_TimerID timerID;	/* The ID to which the interval beginning
                                 * with this marker should be attributed */
    void * marker;
    //cudaEvent_t marker; 		/* The driver event for this marker */
    struct pb_async_time_marker_list *next;
};

struct pb_SubTimer {
    char *label;
    struct pb_Timer timer;
    struct pb_SubTimer *next;
};

struct pb_SubTimerList {
    struct pb_SubTimer *current;
    struct pb_SubTimer *subtimer_list;
};

/* A set of timers for recording execution times. */
struct pb_TimerSet {
    enum pb_TimerID current;
    struct pb_async_time_marker_list* async_markers;
    pb_Timestamp async_begin;
    pb_Timestamp wall_begin;
    struct pb_Timer timers[pb_TimerID_LAST];
    struct pb_SubTimerList *sub_timer_list[pb_TimerID_LAST];
};

/* Reset all timers in the set. */
void
pb_InitializeTimerSet(struct pb_TimerSet *timers);

void
pb_AddSubTimer(struct pb_TimerSet *timers, char *label, enum pb_TimerID pb_Category);

/* Select which timer the next interval of time should be accounted
 * to. The selected timer is started and other timers are stopped.
 * Using pb_TimerID_NONE stops all timers. */
void
pb_SwitchToTimer(struct pb_TimerSet *timers, enum pb_TimerID timer);

void
pb_SwitchToSubTimer(struct pb_TimerSet *timers, char *label, enum pb_TimerID category);

/* Print timer values to standard output. */
void
pb_PrintTimerSet(struct pb_TimerSet *timers, double *time);

/* Release timer resources */
void
pb_DestroyTimerSet(struct pb_TimerSet * timers);

void
pb_SetOpenCL(void *clContextPtr, void *clCommandQueuePtr);

#ifdef __cplusplus
}
#endif

#endif

// end utils.h


// utils.c

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#if _POSIX_VERSION >= 200112L
# include <sys/time.h>
#endif

/* Free an array of owned strings. */
static void
free_string_array(char **string_array)
{
    char **p;

    if (!string_array) return;
    for (p = string_array; *p; p++) free(*p);
    free(string_array);
}

/* Parse a comma-delimited list of strings into an
 * array of strings. */
static char **
read_string_array(char *in)
{
    char **ret;
    int i;
    int count;			/* Number of items in the input */
    char *substring;		/* Current substring within 'in' */

    /* Count the number of items in the string */
    count = 1;
    for (i = 0; in[i]; i++) if (in[i] == ',') count++;

    /* Allocate storage */
    ret = (char **)malloc((count + 1) * sizeof(char *));

    /* Create copies of the strings from the list */
    substring = in;
    for (i = 0; i < count; i++) {
        char *substring_end;
        int substring_length;

        /* Find length of substring */
        for (substring_end = substring;
             (*substring_end != ',') && (*substring_end != 0);
             substring_end++);

        substring_length = substring_end - substring;

        /* Allocate memory and copy the substring */
        ret[i] = (char *)malloc(substring_length + 1);
        memcpy(ret[i], substring, substring_length);
        ret[i][substring_length] = 0;

        /* go to next substring */
        substring = substring_end + 1;
    }
    ret[i] = NULL;		/* Write the sentinel value */

    return ret;
}

struct argparse {
    int argc;			/* Number of arguments.  Mutable. */
    char **argv;			/* Argument values.  Immutable. */

    int argn;			/* Current argument number. */
    char **argv_get;		/* Argument value being read. */
    char **argv_put;		/* Argument value being written.
				 * argv_put <= argv_get. */
};

static void
initialize_argparse(struct argparse *ap, int argc, char **argv)
{
    ap->argc = argc;
    ap->argn = 0;
    ap->argv_get = ap->argv_put = ap->argv = argv;
}

static void
finalize_argparse(struct argparse *ap)
{
    /* Move the remaining arguments */
    for(; ap->argn < ap->argc; ap->argn++)
        *ap->argv_put++ = *ap->argv_get++;
}

/* Delete the current argument. */
static void
delete_argument(struct argparse *ap)
{
    if (ap->argn >= ap->argc) {
        fprintf(stderr, "delete_argument\n");
    }
    ap->argc--;
    ap->argv_get++;
}

/* Go to the next argument.  Also, move the current argument to its
 * final location in argv. */
static void
next_argument(struct argparse *ap)
{
    if (ap->argn >= ap->argc) {
        fprintf(stderr, "next_argument\n");
    }
    /* Move argument to its new location. */
    *ap->argv_put++ = *ap->argv_get++;
    ap->argn++;
}

static int
is_end_of_arguments(struct argparse *ap)
{
    return ap->argn == ap->argc;
}

static char *
get_argument(struct argparse *ap)
{
    return *ap->argv_get;
}

static char *
consume_argument(struct argparse *ap)
{
    char *ret = get_argument(ap);
    delete_argument(ap);
    return ret;
}

struct pb_Parameters *
pb_ReadParameters(int *_argc, char **argv)
{
    char *err_message;
    struct argparse ap;
    struct pb_Parameters *ret =
            (struct pb_Parameters *)malloc(sizeof(struct pb_Parameters));

    /* Initialize the parameters structure */
    ret->outFile = NULL;
    ret->inpFiles = (char **)malloc(sizeof(char *));
    ret->inpFiles[0] = NULL;

    /* Each argument */
    initialize_argparse(&ap, *_argc, argv);
    while(!is_end_of_arguments(&ap)) {
        char *arg = get_argument(&ap);

        /* Single-character flag */
        if ((arg[0] == '-') && (arg[1] != 0) && (arg[2] == 0)) {
            delete_argument(&ap);	/* This argument is consumed here */

            switch(arg[1]) {
                case 'o':			/* Output file name */
                    if (is_end_of_arguments(&ap))
                    {
                        err_message = "Expecting file name after '-o'\n";
                        goto error;
                    }
                    free(ret->outFile);
                    ret->outFile = strdup(consume_argument(&ap));
                    break;
                case 'i':			/* Input file name */
                    if (is_end_of_arguments(&ap))
                    {
                        err_message = "Expecting file name after '-i'\n";
                        goto error;
                    }
                    ret->inpFiles = read_string_array(consume_argument(&ap));
                    break;
                case '-':			/* End of options */
                    goto end_of_options;
                default:
                    err_message = "Unexpected command-line parameter\n";
                    goto error;
            }
        }
        else {
            /* Other parameters are ignored */
            next_argument(&ap);
        }
    } /* end for each argument */

    end_of_options:
    *_argc = ap.argc;		/* Save the modified argc value */
    finalize_argparse(&ap);

    return ret;

    error:
    fputs(err_message, stderr);
    pb_FreeParameters(ret);
    return NULL;
}

void
pb_FreeParameters(struct pb_Parameters *p)
{
    char **cpp;

    free(p->outFile);
    free_string_array(p->inpFiles);
    free(p);
}

int
pb_Parameters_CountInputs(struct pb_Parameters *p)
{
    int n;

    for (n = 0; p->inpFiles[n]; n++);
    return n;
}

/*****************************************************************************/
/* Timer routines */

static void
accumulate_time(pb_Timestamp *accum,
                pb_Timestamp start,
                pb_Timestamp end)
{
#if _POSIX_VERSION >= 200112L
    *accum += end - start;
#else
# error "Timestamps not implemented for this system"
#endif
}

#if _POSIX_VERSION >= 200112L
static pb_Timestamp get_time()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (pb_Timestamp) (tv.tv_sec * 1000000LL + tv.tv_usec);
}
#else
# error "no supported time libraries are available on this platform"
#endif

void
pb_ResetTimer(struct pb_Timer *timer)
{
    timer->state = pb_Timer_STOPPED;

#if _POSIX_VERSION >= 200112L
    timer->elapsed = 0;
#else
# error "pb_ResetTimer: not implemented for this system"
#endif
}

void
pb_StartTimer(struct pb_Timer *timer)
{
    if (timer->state != pb_Timer_STOPPED) {
        fputs("Ignoring attempt to start a running timer\n", stderr);
        return;
    }

    timer->state = pb_Timer_RUNNING;

#if _POSIX_VERSION >= 200112L
    {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        timer->init = tv.tv_sec * 1000000LL + tv.tv_usec;
    }
#else
# error "pb_StartTimer: not implemented for this system"
#endif
}

void
pb_StartTimerAndSubTimer(struct pb_Timer *timer, struct pb_Timer *subtimer)
{
    unsigned int numNotStopped = 0x3; // 11
    if (timer->state != pb_Timer_STOPPED) {
        fputs("Warning: Timer was not stopped\n", stderr);
        numNotStopped &= 0x1; // Zero out 2^1
    }
    if (subtimer->state != pb_Timer_STOPPED) {
        fputs("Warning: Subtimer was not stopped\n", stderr);
        numNotStopped &= 0x2; // Zero out 2^0
    }
    if (numNotStopped == 0x0) {
        fputs("Ignoring attempt to start running timer and subtimer\n", stderr);
        return;
    }

    timer->state = pb_Timer_RUNNING;
    subtimer->state = pb_Timer_RUNNING;

#if _POSIX_VERSION >= 200112L
    {
        struct timeval tv;
        gettimeofday(&tv, NULL);

        if (numNotStopped & 0x2) {
            timer->init = tv.tv_sec * 1000000LL + tv.tv_usec;
        }

        if (numNotStopped & 0x1) {
            subtimer->init = tv.tv_sec * 1000000LL + tv.tv_usec;
        }
    }
#else
# error "pb_StartTimer: not implemented for this system"
#endif

}

void
pb_StopTimer(struct pb_Timer *timer)
{

    pb_Timestamp fini;

    if (timer->state != pb_Timer_RUNNING) {
        fputs("Ignoring attempt to stop a stopped timer\n", stderr);
        return;
    }

    timer->state = pb_Timer_STOPPED;

#if _POSIX_VERSION >= 200112L
    {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        fini = tv.tv_sec * 1000000LL + tv.tv_usec;
    }
#else
# error "pb_StopTimer: not implemented for this system"
#endif

    accumulate_time(&timer->elapsed, timer->init, fini);
    timer->init = fini;

}

void pb_StopTimerAndSubTimer(struct pb_Timer *timer, struct pb_Timer *subtimer) {

    pb_Timestamp fini;

    unsigned int numNotRunning = 0x3; // 0b11
    if (timer->state != pb_Timer_RUNNING) {
        fputs("Warning: Timer was not running\n", stderr);
        numNotRunning &= 0x1; // Zero out 2^1
    }
    if (subtimer->state != pb_Timer_RUNNING) {
        fputs("Warning: Subtimer was not running\n", stderr);
        numNotRunning &= 0x2; // Zero out 2^0
    }
    if (numNotRunning == 0x0) {
        fputs("Ignoring attempt to stop stopped timer and subtimer\n", stderr);
        return;
    }


    timer->state = pb_Timer_STOPPED;
    subtimer->state = pb_Timer_STOPPED;

#if _POSIX_VERSION >= 200112L
    {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        fini = tv.tv_sec * 1000000LL + tv.tv_usec;
    }
#else
# error "pb_StopTimer: not implemented for this system"
#endif

    if (numNotRunning & 0x2) {
        accumulate_time(&timer->elapsed, timer->init, fini);
        timer->init = fini;
    }

    if (numNotRunning & 0x1) {
        accumulate_time(&subtimer->elapsed, subtimer->init, fini);
        subtimer->init = fini;
    }

}

/* Get the elapsed time in seconds. */
double
pb_GetElapsedTime(struct pb_Timer *timer)
{
    double ret;

    if (timer->state != pb_Timer_STOPPED) {
        fputs("Elapsed time from a running timer is inaccurate\n", stderr);
    }

#if _POSIX_VERSION >= 200112L
    ret = timer->elapsed / 1e6;
#else
# error "pb_GetElapsedTime: not implemented for this system"
#endif
    return ret;
}

void
pb_InitializeTimerSet(struct pb_TimerSet *timers)
{
    int n;

    timers->wall_begin = get_time();

    timers->current = pb_TimerID_NONE;

    timers->async_markers = NULL;


    for (n = 0; n < pb_TimerID_LAST; n++) {
        pb_ResetTimer(&timers->timers[n]);
        timers->sub_timer_list[n] = NULL; // free first?
    }
}

void
pb_AddSubTimer(struct pb_TimerSet *timers, char *label, enum pb_TimerID pb_Category) {

    struct pb_SubTimer *subtimer = (struct pb_SubTimer *) malloc
            (sizeof(struct pb_SubTimer));

    int len = strlen(label);

    subtimer->label = (char *) malloc (sizeof(char)*(len+1));
    sprintf(subtimer->label, "%s\0", label);

    pb_ResetTimer(&subtimer->timer);
    subtimer->next = NULL;

    struct pb_SubTimerList *subtimerlist = timers->sub_timer_list[pb_Category];
    if (subtimerlist == NULL) {
        subtimerlist = (struct pb_SubTimerList *) malloc
                (sizeof(struct pb_SubTimerList));
        subtimerlist->subtimer_list = subtimer;
        timers->sub_timer_list[pb_Category] = subtimerlist;
    } else {
        // Append to list
        struct pb_SubTimer *element = subtimerlist->subtimer_list;
        while (element->next != NULL) {
            element = element->next;
        }
        element->next = subtimer;
    }

}

void
pb_SwitchToSubTimer(struct pb_TimerSet *timers, char *label, enum pb_TimerID category)
{

// switchToSub( NULL, NONE
// switchToSub( NULL, some
// switchToSub( some, some
// switchToSub( some, NONE -- tries to find "some" in NONE's sublist, which won't be printed

    struct pb_Timer *topLevelToStop = NULL;
    if (timers->current != category && timers->current != pb_TimerID_NONE) {
        // Switching to subtimer in a different category needs to stop the top-level current, different categoried timer.
        // NONE shouldn't have a timer associated with it, so exclude from branch
        topLevelToStop = &timers->timers[timers->current];
    }

    struct pb_SubTimerList *subtimerlist = timers->sub_timer_list[timers->current];
    struct pb_SubTimer *curr = (subtimerlist == NULL) ? NULL : subtimerlist->current;

    if (timers->current != pb_TimerID_NONE) {
        if (curr != NULL && topLevelToStop != NULL) {
            pb_StopTimerAndSubTimer(topLevelToStop, &curr->timer);
        } else if (curr != NULL) {
            pb_StopTimer(&curr->timer);
        } else {
            pb_StopTimer(topLevelToStop);
        }
    }

    subtimerlist = timers->sub_timer_list[category];
    struct pb_SubTimer *subtimer = NULL;

    if (label != NULL) {
        subtimer = subtimerlist->subtimer_list;
        while (subtimer != NULL) {
            if (strcmp(subtimer->label, label) == 0) {
                break;
            } else {
                subtimer = subtimer->next;
            }
        }
    }

    if (category != pb_TimerID_NONE) {

        if (subtimerlist != NULL) {
            subtimerlist->current = subtimer;
        }

        if (category != timers->current && subtimer != NULL) {
            pb_StartTimerAndSubTimer(&timers->timers[category], &subtimer->timer);
        } else if (subtimer != NULL) {
            // Same category, different non-NULL subtimer
            pb_StartTimer(&subtimer->timer);
        } else{
            // Different category, but no subtimer (not found or specified as NULL) -- unprefered way of setting topLevel timer
            pb_StartTimer(&timers->timers[category]);
        }
    }

    timers->current = category;

}

void
pb_SwitchToTimer(struct pb_TimerSet *timers, enum pb_TimerID timer)
{
    /* Stop the currently running timer */
    if (timers->current != pb_TimerID_NONE) {
        struct pb_SubTimer *currSubTimer = NULL;
        struct pb_SubTimerList *subtimerlist = timers->sub_timer_list[timers->current];

        if ( subtimerlist != NULL) {
            currSubTimer = timers->sub_timer_list[timers->current]->current;
        }
        if ( currSubTimer!= NULL) {
            pb_StopTimerAndSubTimer(&timers->timers[timers->current], &currSubTimer->timer);
        } else {
            pb_StopTimer(&timers->timers[timers->current]);
        }

    }

    timers->current = timer;

    if (timer != pb_TimerID_NONE) {
        pb_StartTimer(&timers->timers[timer]);
    }
}

void
pb_PrintTimerSet(struct pb_TimerSet *timers, double *time)
{

    pb_Timestamp wall_end = get_time();

    struct pb_Timer *t = timers->timers;
    struct pb_SubTimer* sub = NULL;

    int maxSubLength;

    const char *categories[] = {
            "IO", "Kernel", "Copy", "Driver", "Copy Async", "Compute"
    };

    const int maxCategoryLength = 10;

    int i;
    for(i = 1; i < pb_TimerID_LAST-1; ++i) { // exclude NONE and OVRELAP from this format
        if(pb_GetElapsedTime(&t[i]) != 0) {

            // Print Category Timer
            printf("%-*s: %f\n", maxCategoryLength, categories[i-1], pb_GetElapsedTime(&t[i]));

            if (timers->sub_timer_list[i] != NULL) {
                sub = timers->sub_timer_list[i]->subtimer_list;
                maxSubLength = 0;
                while (sub != NULL) {
                    // Find longest SubTimer label
                    if (strlen(sub->label) > maxSubLength) {
                        maxSubLength = strlen(sub->label);
                    }
                    sub = sub->next;
                }

                // Fit to Categories
                if (maxSubLength <= maxCategoryLength) {
                    maxSubLength = maxCategoryLength;
                }

                sub = timers->sub_timer_list[i]->subtimer_list;

                // Print SubTimers
                while (sub != NULL) {
                    printf(" -%-*s: %f\n", maxSubLength, sub->label, pb_GetElapsedTime(&sub->timer));
                    sub = sub->next;
                }
            }
        }
    }

    if(pb_GetElapsedTime(&t[pb_TimerID_OVERLAP]) != 0)
        printf("CPU/Kernel Overlap: %f\n", pb_GetElapsedTime(&t[pb_TimerID_OVERLAP]));

    float walltime = (wall_end - timers->wall_begin)/ 1e6;
    printf("Timer Wall Time: %f\n", walltime);

    (*time) = walltime;

}

//void
//pb_PrintTimerSet(struct pb_TimerSet *timers)
//{
//
//    pb_Timestamp wall_end = get_time();
//
//    struct pb_Timer *t = timers->timers;
//    struct pb_SubTimer* sub = NULL;
//
//    int maxSubLength;
//
//    const char *categories[] = {
//            "IO", "Kernel", "Copy", "Driver", "Copy Async", "Compute"
//    };
//
//    const int maxCategoryLength = 10;
//
//    int i;
//    for(i = 1; i < pb_TimerID_LAST-1; ++i) { // exclude NONE and OVRELAP from this format
//        if(pb_GetElapsedTime(&t[i]) != 0) {
//
//            // Print Category Timer
//            printf("%-*s: %f\n", maxCategoryLength, categories[i-1], pb_GetElapsedTime(&t[i]));
//
//            if (timers->sub_timer_list[i] != NULL) {
//                sub = timers->sub_timer_list[i]->subtimer_list;
//                maxSubLength = 0;
//                while (sub != NULL) {
//                    // Find longest SubTimer label
//                    if (strlen(sub->label) > maxSubLength) {
//                        maxSubLength = strlen(sub->label);
//                    }
//                    sub = sub->next;
//                }
//
//                // Fit to Categories
//                if (maxSubLength <= maxCategoryLength) {
//                    maxSubLength = maxCategoryLength;
//                }
//
//                sub = timers->sub_timer_list[i]->subtimer_list;
//
//                // Print SubTimers
//                while (sub != NULL) {
//                    printf(" -%-*s: %f\n", maxSubLength, sub->label, pb_GetElapsedTime(&sub->timer));
//                    sub = sub->next;
//                }
//            }
//        }
//    }
//
//    if(pb_GetElapsedTime(&t[pb_TimerID_OVERLAP]) != 0)
//        printf("CPU/Kernel Overlap: %f\n", pb_GetElapsedTime(&t[pb_TimerID_OVERLAP]));
//
//    float walltime = (wall_end - timers->wall_begin)/ 1e6;
//    printf("Timer Wall Time: %f\n", walltime);
//
//}

void pb_DestroyTimerSet(struct pb_TimerSet * timers)
{
    /* clean up all of the async event markers */
    struct pb_async_time_marker_list ** event = &(timers->async_markers);
    while( *event != NULL) {
        struct pb_async_time_marker_list ** next = &((*event)->next);
        free(*event);
        (*event) = NULL;
        event = next;
    }

    int i = 0;
    for(i = 0; i < pb_TimerID_LAST; ++i) {
        if (timers->sub_timer_list[i] != NULL) {
            struct pb_SubTimer *subtimer = timers->sub_timer_list[i]->subtimer_list;
            struct pb_SubTimer *prev = NULL;
            while (subtimer != NULL) {
                free(subtimer->label);
                prev = subtimer;
                subtimer = subtimer->next;
                free(prev);
            }
            free(timers->sub_timer_list[i]);
        }
    }
}

// end utils.c


// dump.h

/***************************************************************************
 *
 *            (C) Copyright 2010 The Board of Trustees of the
 *                        University of Illinois
 *                         All Rights Reserved
 *
 ***************************************************************************/

void dump_histo_img(unsigned char* histo, unsigned int height, unsigned int width, const char *filename);

//end dump.h


// dump.c

/***************************************************************************
 *
 *            (C) Copyright 2010 The Board of Trustees of the
 *                        University of Illinois
 *                         All Rights Reserved
 *
 ***************************************************************************/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


// This function takes an HSV value and converts it to BMP.
// We use this function to generate colored images with
// Smooth spectrum traversal for the input and output images.
RGB HSVtoRGB( float h, float s, float v )
{
    int i;
    float f, p, q, t;
    float r, g, b;
    RGB value={0,0,0};

    if( s == 0 ) {
        r = g = b = v;
        return value;
    }
    h /= 60;
    i = floor( h );
    f = h - i;
    p = v * ( 1 - s );
    q = v * ( 1 - s * f );
    t = v * ( 1 - s * ( 1 - f ) );
    switch( i ) {
        case 0:
            r = v; g = t; b = p;
            break;
        case 1:
            r = q; g = v; b = p;
            break;
        case 2:
            r = p; g = v; b = t;
            break;
        case 3:
            r = p; g = q; b = v;
            break;
        case 4:
            r = t; g = p; b = v;
            break;
        default:
            r = v; g = p; b = q;
            break;
    }

    unsigned int temp = r*255;
    value.R = temp;
    temp = g*255;
    value.G = temp;
    temp = b*255;
    value.B = temp;

    return value;
}

void dump_histo_img(unsigned char* histo, unsigned int height, unsigned int width, const char *filename)
{
    RGB* pixel_map = (RGB*) malloc (height*width*sizeof(RGB));

    size_t y, x;
    for (y = 0; y < height; ++y)
    {
        for (x = 0; x < width; ++x)
        {
            unsigned char value = histo[y * width + x];

            if (value == 0){
                pixel_map[y*width+x].R = 0;
                pixel_map[y*width+x].G = 0;
                pixel_map[y*width+x].B = 0;
            } else {
                pixel_map[y*width+x] = HSVtoRGB(0.0,1.0,cbrt(1+ 63.0*((float)value)/((float)UINT8_MAX))/4);
            }
        }
    }
    create_bmp(pixel_map, height, width, filename);
    free(pixel_map);
}

// end dump.c


// main.c

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "common.h"

typedef struct {

} hist;

int sequential(
        int argc,
        char* argv[],
        struct pb_Parameters *parameters,
        unsigned char **histo,
        unsigned int *size,
        double *time
) {
    printf("\nSequential execution:\n");

    struct pb_TimerSet timers;

    if(!parameters->inpFiles[0]){
        fputs("Input file expected\n", stderr);
        return -1;
    }

    int numIterations;
    if (argc >= 2){
        numIterations = atoi(argv[1]);
    } else {
        fputs("Expected at least one command line argument\n", stderr);
        return -1;
    }

    pb_InitializeTimerSet(&timers);

    char *inputStr = "Input";
    char *outputStr = "Output";

    pb_AddSubTimer(&timers, inputStr, pb_TimerID_IO);
    pb_AddSubTimer(&timers, outputStr, pb_TimerID_IO);

    pb_SwitchToSubTimer(&timers, inputStr, pb_TimerID_IO);

    unsigned int img_width, img_height;
    unsigned int histo_width, histo_height;

    FILE* f = fopen(parameters->inpFiles[0],"rb");
    int result = 0;

    result += fread(&img_width,    sizeof(unsigned int), 1, f);
    result += fread(&img_height,   sizeof(unsigned int), 1, f);
    result += fread(&histo_width,  sizeof(unsigned int), 1, f);
    result += fread(&histo_height, sizeof(unsigned int), 1, f);

    if (result != 4){
        fputs("Error reading input and output dimensions from file\n", stderr);
        return -1;
    }

    unsigned int histo_size = histo_width*histo_height;
    unsigned int img_size = img_width*img_height;
    (*size) = histo_size;

    unsigned int* img = (unsigned int*) malloc (img_size * sizeof(unsigned int));
    (*histo) = (unsigned char*) calloc (histo_size, sizeof(unsigned char));

    pb_SwitchToSubTimer(&timers, "Input", pb_TimerID_IO);

    result = fread(img, sizeof(unsigned int), img_width*img_height, f);

    fclose(f);

    if (result != img_width*img_height){
        fputs("Error reading input array from file\n", stderr);
        return -1;
    }

    pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);


    int iter;
    for (iter = 0; iter < numIterations; iter++){
        memset((*histo),0,histo_height*histo_width*sizeof(unsigned char));
        unsigned int i;
        for (i = 0; i < img_size; ++i) {
            const unsigned int value = img[i];
            if ((*histo)[value] < UINT8_MAX) {
                ++(*histo)[value];
            }
        }
    }

    pb_SwitchToSubTimer(&timers, outputStr, pb_TimerID_IO);

    if (parameters->outFile) {
        dump_histo_img((*histo), histo_height, histo_width, parameters->outFile);
    }

    pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);

    free(img);

    pb_SwitchToTimer(&timers, pb_TimerID_NONE);

    printf("\n");
    pb_PrintTimerSet(&timers, time);

    return 0;
}

int parallel(
        int argc,
        char* argv[],
        struct pb_Parameters *parameters,
        unsigned char **histo,
        unsigned int *size,
        double *time
) {
    struct pb_TimerSet timers;

    printf("\nParallel execution:\n");

    if(!parameters->inpFiles[0]){
        fputs("Input file expected\n", stderr);
        return -1;
    }

    int numIterations;
    if (argc >= 2){
        numIterations = atoi(argv[1]);
    } else {
        fputs("Expected at least one command line argument\n", stderr);
        return -1;
    }

    pb_InitializeTimerSet(&timers);

    char *inputStr = "Input";
    char *outputStr = "Output";

    pb_AddSubTimer(&timers, inputStr, pb_TimerID_IO);
    pb_AddSubTimer(&timers, outputStr, pb_TimerID_IO);

    pb_SwitchToSubTimer(&timers, inputStr, pb_TimerID_IO);

    unsigned int img_width, img_height;
    unsigned int histo_width, histo_height;

    FILE* f = fopen(parameters->inpFiles[0],"rb");
    int result = 0;

    result += fread(&img_width,    sizeof(unsigned int), 1, f);
    result += fread(&img_height,   sizeof(unsigned int), 1, f);
    result += fread(&histo_width,  sizeof(unsigned int), 1, f);
    result += fread(&histo_height, sizeof(unsigned int), 1, f);

    if (result != 4){
        fputs("Error reading input and output dimensions from file\n", stderr);
        return -1;
    }

    unsigned int histo_size = histo_height * histo_width;
    unsigned int img_size = img_width * img_height;
    *size = histo_size;

    unsigned int* img = (unsigned int*) malloc (img_width*img_height*sizeof(unsigned int));
    (*histo) = (unsigned char*) calloc (histo_width*histo_height, sizeof(unsigned char));

    pb_SwitchToSubTimer(&timers, "Input", pb_TimerID_IO);

    result = fread(img, sizeof(unsigned int), img_width*img_height, f);

    fclose(f);

    if (result != img_width*img_height){
        fputs("Error reading input array from file\n", stderr);
        return -1;
    }

    pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);

    //int iter;

#pragma omp parallel
    {
        unsigned char *inner_histo = (unsigned char *) calloc(histo_width * histo_height, sizeof(unsigned char));

        double num_threads = omp_get_num_threads();
        int thread_id = omp_get_thread_num();

        int chunk_size = ceil(numIterations/num_threads);
        int iter_start = thread_id * chunk_size;
        int iter_end = (thread_id + 1) * chunk_size;

        for (int iter = iter_start; iter < iter_end; iter++) {

                memset(inner_histo, 0, histo_height * histo_width * sizeof(unsigned char));
                unsigned int i;
                for (i = 0; i < img_width * img_height; ++i) {

#pragma omp task
                        {
                        const unsigned int value = img[i];
                        if (inner_histo[value] < UINT8_MAX) {
                            ++inner_histo[value];
                        }
                    }
                }

                for (i = 0; i < histo_size; i++) {
                    (*histo)[i] = inner_histo[i];
                }
        }

        free(inner_histo);
    }

    pb_SwitchToSubTimer(&timers, outputStr, pb_TimerID_IO);

    if (parameters->outFile) {
        dump_histo_img((*histo), histo_height, histo_width, parameters->outFile);
    }

    pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);

    free(img);

    pb_SwitchToTimer(&timers, pb_TimerID_NONE);

    printf("\n");
    pb_PrintTimerSet(&timers, time);

    return 0;
}

//int parallel(
//        int argc,
//        char* argv[],
//        struct pb_Parameters *parameters,
//        unsigned char **histo,
//        unsigned int *size,
//        double *time
//) {
//    struct pb_TimerSet timers;
//
//    printf("\nParallel execution:\n");
//
//    if(!parameters->inpFiles[0]){
//        fputs("Input file expected\n", stderr);
//        return -1;
//    }
//
//    int numIterations;
//    if (argc >= 2){
//        numIterations = atoi(argv[1]);
//    } else {
//        fputs("Expected at least one command line argument\n", stderr);
//        return -1;
//    }
//
//    pb_InitializeTimerSet(&timers);
//
//    char *inputStr = "Input";
//    char *outputStr = "Output";
//
//    pb_AddSubTimer(&timers, inputStr, pb_TimerID_IO);
//    pb_AddSubTimer(&timers, outputStr, pb_TimerID_IO);
//
//    pb_SwitchToSubTimer(&timers, inputStr, pb_TimerID_IO);
//
//    unsigned int img_width, img_height;
//    unsigned int histo_width, histo_height;
//
//    FILE* f = fopen(parameters->inpFiles[0],"rb");
//    int result = 0;
//
//    result += fread(&img_width,    sizeof(unsigned int), 1, f);
//    result += fread(&img_height,   sizeof(unsigned int), 1, f);
//    result += fread(&histo_width,  sizeof(unsigned int), 1, f);
//    result += fread(&histo_height, sizeof(unsigned int), 1, f);
//
//    if (result != 4){
//        fputs("Error reading input and output dimensions from file\n", stderr);
//        return -1;
//    }
//
//    unsigned int histo_size = histo_height * histo_width;
//    unsigned int img_size = img_width*img_height;
//    (*size) = histo_size;
//
//    unsigned int* img = (unsigned int*) malloc (img_width*img_height*sizeof(unsigned int));
//    (*histo) = (unsigned char*) calloc (histo_width*histo_height, sizeof(unsigned char));
//
//    pb_SwitchToSubTimer(&timers, "Input", pb_TimerID_IO);
//
//    result = fread(img, sizeof(unsigned int), img_width*img_height, f);
//
//    fclose(f);
//
//    if (result != img_width*img_height){
//        fputs("Error reading input array from file\n", stderr);
//        return -1;
//    }
//
//    pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);
//
//    //int iter;
//
//    for (int iter = 0; iter < numIterations; iter++) {
//
//        memset((*histo), 0, histo_size * sizeof(unsigned char));
//
//        unsigned char *hista = (unsigned char*) malloc (histo_size *4 * sizeof(unsigned char));
//
//        //memset(hista, 0, histo_size * 4);
//
//    int thread_id, offset, iter_start, iter_end, num_thread, chunk_size;
//#pragma omp parallel \
//    private (thread_id, offset, iter_start, iter_end) \
//    shared (chunk_size, hista, num_thread)
//        {
//            num_thread = omp_get_num_threads();
//            thread_id = omp_get_thread_num();
//            offset = thread_id * histo_size;
//
//            memset(hista + offset, 0 , histo_size);
//
//            chunk_size = ceil(img_size/(double)num_thread);
//            iter_start = thread_id * chunk_size;
//            iter_end = (thread_id + 1) * chunk_size;
//            if (iter_end > img_size) iter_end = img_size;
////#pragma omp critical
//        //{
//            for (int i = iter_start; i < iter_end; ++i) {
//                const unsigned int value = img[i];
//                if (hista[offset + value] < UINT8_MAX) {
//                    ++hista[offset + value];
//                }
//            }
//        //}
//#pragma omp barrier
//        }
//
//
//#pragma omp for
//        for(int i=0; i<histo_size; i++) {
//            for(int t=0; t<num_thread; t++) {
//                if ((*histo)[i] < UINT8_MAX - hista[histo_size*t + i]) {
//                    (*histo)[i] += hista[histo_size*t + i];
//                }
//                else {
//                    (*histo)[i] = UINT8_MAX;
//                }
//            }
//        }
//
//        free(hista);
//    }
//
//    pb_SwitchToSubTimer(&timers, outputStr, pb_TimerID_IO);
//
//    if (parameters->outFile) {
//        dump_histo_img((*histo), histo_height, histo_width, parameters->outFile);
//    }
//
//    pb_SwitchToTimer(&timers, pb_TimerID_COMPUTE);
//
//    free(img);
//
//    pb_SwitchToTimer(&timers, pb_TimerID_NONE);
//
//    printf("\n");
//    pb_PrintTimerSet(&timers, time);
//
//    return 0;
//}


int main(int argc, char* argv[]) {
    double sequential_time, parallel_time;
    int err;

    unsigned char* sequential_histo, *parallel_histo;
    unsigned int sequential_histo_size, parallel_histo_size;

    struct pb_Parameters *parameters;
    parameters = pb_ReadParameters(&argc, argv);
    if (!parameters)
        return -1;

    err = parallel(argc, argv, parameters, &parallel_histo, &parallel_histo_size, &parallel_time);
    if (err) { return err; }

    err = sequential(argc, argv, parameters, &sequential_histo, &sequential_histo_size, &sequential_time);
    if (err) { return err; }

    finish_2(sequential_histo, parallel_histo, sequential_histo_size, parallel_histo_size, sequential_time, parallel_time);

    pb_FreeParameters(parameters);

    free(sequential_histo);
    free(parallel_histo);
}

// end main.c