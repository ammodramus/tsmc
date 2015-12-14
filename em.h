typedef struct 
{
    Data * dat;
    Hmm * hmm;
    int numSeqs;
    SeqType seqtype;
    double *** forward;
    double *** backward;
    double ** firstPosterior;
    double ** expect;
    double ** normConst;
} Em;

void Em_init(Em * em, Data * dat, Hmm * hmm);
void Em_free(Em * em);
void Em_get_forward(Em * em);
