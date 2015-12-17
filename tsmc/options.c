#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "definitions.h"
#include "options.h"


static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
    {"iterations", required_argument, 0, 'i'},
    {"filename", required_argument, 0, 'f'},
    {"popsizes", required_argument, 0, 'p'},
	{0, 0, 0, 0}
};

char usage[] = "\ntsmc [OPTIONS]\n\
    \n\
	\n";

void Options_print_help_statement()
{
	printf("%s", &usage[0]);
	return;
}

void Options_set_defaults(Options * opt)
{
    strcpy(&(opt->filename[0]), "-");
    // 2*4+8*2+1*4
    int i;
    opt->numFreeLambdas = 10;
    int counts[] = {4,4,2,2,2,2,2,2,2,2,4};
    int sum = 0;
    for(i = 0; i < opt->numFreeLambdas+1; i++)
    {
        opt->lambdaCounts[i] = counts[i];
        sum += counts[i];
    }
    opt->n = sum;

    opt->numEmIterations = 20;

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

/*
 * Algorithm for parsing options
 *
 * initialize
 * -----
 *  expectNum = True
 *  times = 1
 *  numIdx = 0
 *  num = empty string
 * 
 * iterate
 * ------
 *  - get character ch
 *  - end if ch is '\0'
 *  - if expectNum:
 *       is not numeric:
 *         error
 *       else:
 *         num[numIdx] = 
 */     

void Options_get_lambda_params(Options * opt)
{
    char * const filename = opt->filename;
    char ch;
    int i;
    char chnum[20];
    int num, success;
    int strIdx = 0, numIdx = 0;
    int done = 0, good;
    int expectNum = 1, times = 1;
    readState state = count;
    opt->numFreeLambdas = 0;
    while((ch = filename[strIdx++]) != '\0')
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
                num[numIdx] = '\0';
                numIdx = 0;
                success = sscanf("%i", &num);
                if(!success)
                {
                    perror("invalid -p argument");
                }
                if(ch == '+')
                {
                    opt->lambdaCounts[opt->numFreeLambdas] = times;
                    for(i = 0; i < times; i++)
                    {
                        opt->lambdaCounts[opt->numFreeLambdas] = num;
                        opt->numFreeLambdas++;
                    }
                    times = 1;
                    expectNum = 1;
                    continue;
                }
                else if(ch == '*')
                {
                    // times
                    times = num;
                    expectNum = 1;
                }
                else
                {
                    perror("invalid character in -p argument");
                }
            }
        }
    }
    return;
}

void Options_parse_options(int argc, char ** argv, Options * opt)
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
		c = getopt_long(argc, argv, "hi:f:p:", long_options, &optionIndex);
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
            case 'f':
                success = sscanf(optarg, "%s", &(opt->filename));
                break;
            case 'p':
                success = sscanf(optarg, "%s", &(opt->paramString));
                break;
			default:
				Options_print_help_statement();
				exit(1);
				break;
		}
	}

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
    assert(opt->numFreeLambdas > 0);
    for(i = 0; i < opt->numFreeLambdas; i++)
    {
        printf("lambda count %i = %i\n", opt->lambdaCounts[i]);
    }
    return;
}
