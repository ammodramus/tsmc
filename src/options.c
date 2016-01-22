#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "definitions.h"
#include "options.h"

#define DEFAULTPATTERN "4+5*3+4"


static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
    {"iterations", required_argument, 0, 'i'},
    {"num-optimizations", required_argument, 0, 'o'},
    {"no-asex", no_argument, 0, 1},
    {"psmc-intervals", no_argument, 0, 2},
    {"psmc-max-time", required_argument, 0, 3},
    {"dip-trip", no_argument, 0, 4},
	{0, 0, 0, 0}
};

char usage[] = "\nUsage: tsmc [OPTIONS] FILE\n\
Performs SMC inference on unphased triploid genome represented in FILE\n\
\n\
    -h, --help                   Display this help message.\n\
    -p=PATT                      Set free population parameter pattern to patt (like PSMC) [4+5*3+4]\n\
    -i, --iterations=ITER        Set number of EM iterations to iter [20]\n\
    -o, --num-optimizations=NUM  Set number of optimizations in maximization step to num [20]\n\
    --dip-trip                   Enable inference of diploidy-to-triploidy transition phase [disabled]\n\
	\n";

void Options_print_help_statement()
{
	printf("%s", &usage[0]);
	return;
}

void Options_set_defaults(Options * opt)
{
    strcpy(&(opt->filename[0]), "-");
    int i;
    opt->numFreeLambdas = 10;
    strcpy(&(opt->paramString[0]), DEFAULTPATTERN); 

    opt->numEmIterations = 20;
    opt->numOptimizations = 20;
    opt->asexEnabled = 1;
    opt->psmcIntervals = 1;
    opt->maxT = 5.0;
    opt->flagDt = 0;

    return;
}

typedef enum
{
    count = 0,
    length
} readState;

const char validChars[] = {'+','*','0','1','2','3','4','5','6','7','8','9'};
const int numValidChars = 12;

int check_char(char ch)
{
    int i;
    for(i = 0; i < numValidChars; i++)
    {
        if(ch == validChars[i])
        {
            return 1;
        }
    }
    return 0;
}

inline int is_numeric(char ch)
{
    char numCh = ch - '0';
    return (numCh >= 0 && numCh <= 9);
}

void Options_get_lambda_params(Options * opt)
{
    char * const pattern = opt->paramString;
    char ch;
    int i;
    char chnum[20];
    int num, success;
    int strIdx = 0, numIdx = 0;
    int done = 0, good;
    int expectNum = 1, times = 1, haveNum = 0;
    readState state = count;
    opt->numFreeLambdas = 0;
    while((ch = pattern[strIdx++]) != '\0')
    {
        good = check_char(ch);
        if(!good)
        {
            perror("invalid character in -p argument");
        }
        if(expectNum)
        {
            if(!is_numeric(ch))
            {
                perror("non-number where expected number in -p argument");
            }
            chnum[numIdx++] = ch;
            expectNum = 0;
            haveNum = 1;
            continue;
        }
        else
        {
            assert(!expectNum);
            if(is_numeric(ch))
            {
                chnum[numIdx++] = ch;
                continue;
            }
            else
            {
                assert(!is_numeric(ch));
                chnum[numIdx] = '\0';
                numIdx = 0;
                success = sscanf(chnum, "%i", &num);
                if(!success)
                {
                    perror("invalid -p argument");
                }
                if(ch == '+')
                {
                    for(i = 0; i < times; i++)
                    {
                        opt->lambdaCounts[opt->numFreeLambdas] = num;
                        opt->numFreeLambdas++;
                    }
                    times = 1;
                    expectNum = 1;
                    haveNum = 0;
                    continue;
                }
                else if(ch == '*')
                {
                    // times
                    times = num;
                    expectNum = 1;
                    haveNum = 0;
                }
                else
                {
                    perror("invalid character in -p argument");
                }
            }
        }
    }
    if(haveNum) // (haveNum might always be equal to !expectNum)
    {
        chnum[numIdx] = '\0';
        success = sscanf(chnum, "%i", &num);
        if(!success)
        {
            perror("invalid -p argument\n");
        }
        for(i = 0; i < times; i++)
        {
            opt->lambdaCounts[opt->numFreeLambdas] = num;
            opt->numFreeLambdas++;
        }
    }
    else
    {
        perror("invalid -p argument");
    }
    opt->numFreeLambdas--;  // the first lambda is not a free param (always
                            // equal to 1)
    opt->n = 0;
    for(i = 0; i < opt->numFreeLambdas+1; i++)
    {
        opt->n += opt->lambdaCounts[i];
    }
    opt->n--;

    return;
}

void Options_parse_options(Options * opt, int argc, char ** argv)
{
	int optionIndex, success;
	char c;
	// set default values
    Options_set_defaults(opt);

    if(argc == 1)
    {
        Options_print_help_statement();
        exit(0);
    }


	while(1)
	{
		c = getopt_long(argc, argv, "hi:p:o:", long_options, &optionIndex);
		if(c == -1)
			break;
		switch(c)
		{
			case 'h':
				Options_print_help_statement();
                exit(0);
				break;
            case 'i':
                success = sscanf(optarg, "%i", &(opt->numEmIterations));
                break;
            case 'o':
                success = sscanf(optarg, "%i", &(opt->numOptimizations));
                break;
            case 'p':
                success = sscanf(optarg, "%s", &(opt->paramString[0]));
                break;
            case 1: // --no-asex
                opt->asexEnabled = 0;
                break;
            case 2: // --psmc-intervals
                opt->psmcIntervals = 1;
                break;
            case 3: // --psmc-max-time
                success = sscanf(optarg, "%lf", &(opt->maxT));
                if(opt->maxT <= 0.0)
                {
                    perror("--psmc-max-time must be > 0");
                }
                break;
            case 4: // --dip-trip
                opt->flagDt = 1;
                break;
			default:
				Options_print_help_statement();
				exit(1);
				break;
		}
	}
    if(optind == argc)
    {
        fprintf(stderr, "\nNo input file specified\n\n");
        Options_print_help_statement();
        exit(1);
    }
    success = sscanf(argv[optind], "%s", &(opt->filename[0]));

    Options_get_lambda_params(opt);

	return;
}

void Options_print_options(Options * opt)
{
    assert(opt);
    int i;
    printf("filename: %s\n", opt->filename);
    printf("numIterations: %i\n", opt->numEmIterations);
    printf("numFreeLambdas: %i\n", opt->numFreeLambdas);
    printf("n: %i\n", opt->n);
    assert(opt->numFreeLambdas >= 0);
    for(i = 0; i < opt->numFreeLambdas+1; i++)
    {
        printf("lambda count %i = %i\n", i, opt->lambdaCounts[i]);
    }
    printf("asex enabled: %i\n", opt->asexEnabled);
    printf("dip-trip enabled: %i\n", opt->flagDt);
    return;
}
