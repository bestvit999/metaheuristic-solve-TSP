#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <time.h>
#include <math.h>
#include <algorithm>

typedef struct node
{
    int index;
    int x;
    int y;
    bool Rbit; // redundant bit
} node;

typedef std::vector<node> nodes;                    // store all nodes information
typedef std::vector<node> trail;                    // the trail path that ant traveled
typedef std::vector<trail> trails;                  // all trails

/* main function */
void main_pmx();
void main_cx();
void main_ox();

/* GA */
// init
trails construct(nodes g); // construct the trails (population)

// three crossover operations
void pmx(trail &chr1, trail &chr2);
void cx(trail &chr1, trail &chr2);
void ox(trail &chr1, trail &chr2);

// function for PMX (partially mapped crossover)
typedef std::vector<nodes> MR; // mapping relationship table
MR exchange_substring(trail &chr1, trail &chr2);
void update_MR(MR &mr, node gene1,node gene2);
void legalize(MR mr,trail &pre_chr1, trail &pre_chr2);
void isRedundant(trail &chromosome, node gene);
bool checkMR(node gene1, node gene2);
node legalize_gene(MR mr, trail chromosome, node gene);
node LS_legal_gene(nodes relation, trail chromosome);
void showRbit(trail chromesome);

// function for CX (cycle crossover)
void find_cycle(trail &chr1, trail &chr2);
int find_gene(trail chromosome, node gene);
void legalize_cx(trail &chr1, trail &chr2);

//function for OX (order crossover)
void find_subtour(trail &chr1, trail &chr2);
void legalize_order(trail &chr1, trail &chr2);
node LS_legalize_order(trail illegal_chr, trail candidate_chr);

// evaluation + determine
trail local_best(trails chromosomes);
trail pick_better_chromosome(trails pop_p);

// mutate
void mutate(trails &chromosomes);
void TwoOpt(trail &chromosome);
trail TwoOptSwap(trail chromosome, const int& i, const int& k );

/* others */
nodes readtsp(std::string path_to_file);
void outtsp(trail t);
double evaluate_weight(node i, node j);
double evaluate_trail_distance(trail t);
bool is_opt(trail best);
void plot(trail best);
void showTrail(trail best);

#define pop_size 10 // population size = 10
#define limited_iter 1250 // end condiction
 
int main(int argc, char *argv[]){
    if (std::stoi(argv[1]) == 1){
        std::cout << "pmx is using" << std::endl;
        main_pmx();
    }else if (std::stoi(argv[1]) == 2){
        std::cout << "cx is using" << std::endl;
        main_cx();
    }else if (std::stoi(argv[1]) == 3){
        std::cout << "ox is using" << std::endl;
        main_ox();
    }
}

void main_pmx()
{
    srand(time(NULL));

    nodes g = readtsp("eil51.txt");
    trails pop_p = construct(g); // population of parent set
    trails pop_c; // population of child set
    trail best = pop_p[0]; // tracking the best trail

    int iter = 0;
    while(iter < limited_iter){

        // evalaute the local_best in population
        trail lbest = local_best(pop_p);
        if (evaluate_trail_distance(best) > evaluate_trail_distance(lbest)){
            best = lbest;
        }

        // crossover with CR to generate the child population
        for (int i = 0;i < pop_p.size() / 2;i++){

            // pick the best kth individuals in parent set to generate pop_c
            // trail a = pick_better_chromosome(pop_p);
            // trail b = pick_better_chromosome(pop_p);

            trail a = pop_p[rand() % pop_p.size()];
            trail b = pop_p[rand() % pop_p.size()];

            pmx(a,b);

            pop_c.push_back(a);
            pop_c.push_back(b);

        }

        pop_p = pop_c;
        pop_c.clear();  

        // mutate all individual in population (kind of local-search)
        mutate(pop_p);
        
        showTrail(best);
        std::cout << evaluate_trail_distance(best) << std::endl;
        plot(best);


        iter++;
    }
    // std::cout << "no bug" << std::endl;
}

void main_cx(){
    srand(time(NULL));

    nodes g = readtsp("eil51.txt");
    trails pop_p = construct(g); // population of parent set
    trails pop_c; // population of child set
    trail best = pop_p[0]; // tracking the best trail

    int iter = 0;
    while(iter < limited_iter){
        // evalaute the local_best in population
        trail lbest = local_best(pop_p);
        if (evaluate_trail_distance(best) > evaluate_trail_distance(lbest)){
            best = lbest;
        }

        // crossover with CR to generate the child population
        for (int i = 0;i < pop_p.size() / 2;i++){

            // pick the best kth individuals in parent set to generate pop_c
            // trail a = pick_better_chromosome(pop_p);
            // trail b = pick_better_chromosome(pop_p);

            trail a = pop_p[rand() % pop_p.size()];
            trail b = pop_p[rand() % pop_p.size()];

            cx(a,b);

            pop_c.push_back(a);
            pop_c.push_back(b);

        }

        pop_p = pop_c;
        pop_c.clear();  

        // mutate all individual in population (kind of local-search)
        mutate(pop_p);
        
        showTrail(best);
        std::cout << evaluate_trail_distance(best) << std::endl;
        plot(best);


        iter++;
    }
}

void main_ox(){
    srand(time(NULL));

    nodes g = readtsp("eil51.txt");
    trails pop_p = construct(g); // population of parent set
    trails pop_c; // population of child set
    trail best = pop_p[0]; // tracking the best trail

    int iter = 0;
    while(iter < limited_iter){
        // evalaute the local_best in population
        trail lbest = local_best(pop_p);
        if (evaluate_trail_distance(best) > evaluate_trail_distance(lbest)){
            best = lbest;
        }

        // crossover with CR to generate the child population
        for (int i = 0;i < pop_p.size() / 2;i++){

            // pick the best kth individuals in parent set to generate pop_c
            // trail a = pick_better_chromosome(pop_p);
            // trail b = pick_better_chromosome(pop_p);

            trail a = pop_p[rand() % pop_p.size()];
            trail b = pop_p[rand() % pop_p.size()];

            ox(a,b);

            pop_c.push_back(a);
            pop_c.push_back(b);

        }

        pop_p = pop_c;
        pop_c.clear();  

        // mutate all individual in population (kind of local-search)
        mutate(pop_p);
        
        showTrail(best);
        std::cout << evaluate_trail_distance(best) << std::endl;
        plot(best);


        iter++;
    }
}


trails construct(nodes g)
{
// construct the trails (population)
// here we just generate each trail ramdomly

    trails pop;
    trail t;


    for (int i = 0;i < pop_size;i++){
    // generate the number of 'pop_size' trails 
        nodes node_candidate = g; // used to generate one trail
    
        while (!node_candidate.empty()){
            int index = rand() % node_candidate.size();
            t.push_back(node_candidate[index]);
            node_candidate.erase(node_candidate.begin()+index);
        }
        pop.push_back(t);
        t.clear();
    }
    return pop;
}

trail local_best(trails chromosomes){
// select the best chromosome in population
    trail local_best = chromosomes[0];
    for (int i = 1;i < chromosomes.size();i++){
        if (evaluate_trail_distance(local_best) > evaluate_trail_distance(chromosomes[i])){
            local_best = chromosomes[i];
        }
    }
    return local_best;
}

void pmx(trail &chr1, trail &chr2){
    MR mr = exchange_substring(chr1, chr2);

    legalize(mr,chr1,chr2);
}

MR exchange_substring(trail &chr1, trail &chr2){
    int effect_range = 7;
    int index_start = rand() % chr1.size();

    MR mr;
    node tmp;

    for (int i = 0;i < effect_range;i++){
        tmp = chr1[(index_start+i) % chr1.size()];
        chr1[(index_start+i) % chr1.size()] = chr2[(index_start+i) % chr2.size()];
        chr2[(index_start+i) % chr2.size()] = tmp;
        // simultaneously mapping relationship
        update_MR(mr,chr1[(index_start+i) % chr1.size()], chr2[(index_start+i) % chr2.size()]);

        // find the redundant
        isRedundant(chr1,chr1[(index_start+i) % chr1.size()]);
        isRedundant(chr2,chr2[(index_start+i) % chr2.size()]);

        // selected range will not be set to redundant
        chr1[(index_start+i) % chr1.size()].Rbit = 0;
        chr2[(index_start+i) % chr2.size()].Rbit = 0;
    }

    return mr;
}

void update_MR(MR &mr, node gene1,node gene2){
    nodes relation;

    bool stop = 0;
    for (int i = 0;i < mr.size() && !stop;i++){
        for (int j = 0;j < mr[i].size() && !stop;j++){
            if (mr[i][j].index == gene1.index){
                mr[i].push_back(gene2);
                stop = 1;
            }else if (mr[i][j].index == gene2.index){
                mr[i].push_back(gene1);
                stop = 1;
            }
        }
    }

    if (!stop){
        relation.push_back(gene1);
        relation.push_back(gene2);
        mr.push_back(relation);
        relation.clear();
    }
    
}

void isRedundant(trail &chromosome, node gene){
// set the redundant gene's R_bit to 1

    for (int i = 0;i < chromosome.size();i++){
        if (chromosome[i].index == gene.index){
            chromosome[i].Rbit = 1;
        }
    }
}

void legalize(MR mr, trail &pre_chr1, trail &pre_chr2){
// legalize the Rbit == 1 genes with Mapping Relationship table

    for (int i = 0;i < pre_chr1.size();i++){
        if (pre_chr1[i].Rbit == 1){
            pre_chr1[i] = legalize_gene(mr, pre_chr1, pre_chr1[i]);
            pre_chr1[i].Rbit = 0;
        }

        if (pre_chr2[i].Rbit == 1){
            pre_chr2[i] = legalize_gene(mr, pre_chr2, pre_chr2[i]);
            pre_chr2[i].Rbit = 0;
        }
    }
}

node legalize_gene(MR mr,trail chromosome, node gene){
// transform to legal_gene by refer to MR
    for (int i = 0;i < mr.size();i++){
        for (int j = 0;j < mr[i].size();j++){
            // find the relation
            if (mr[i][j].index == gene.index){
                return LS_legal_gene(mr[i],chromosome);
            }
        }
    }
    return gene;
}

node LS_legal_gene(nodes relation, trail chromosome){
// locally find the relationship in the chromosome
    int index = 0;
    for (int i = 0;i < relation.size();i++){
        int count = 0;
        for (int j = 0;j < chromosome.size();j++){
            if (relation[i].index == chromosome[j].index){
                count++;
            }
        }
        if (count == 0){
            index = i;
            return relation[i];
        }
    }
    return relation[index];
}

void cx(trail &chr1, trail &chr2){

    find_cycle(chr1,chr2);
    legalize_cx(chr1,chr2);

}

void find_cycle(trail &chr1, trail &chr2){
    // find the cycle between two chromosome, the rest genes' Rbit will be mark to 1

    for (int i = 0; i < chr1.size();i++){
        chr1[i].Rbit = 1;
        chr2[i].Rbit = 1;
    }

    // find the cycle in chr1
    node start_node = chr1[0];
    chr1[0].Rbit = 0;
    node next = start_node;

    int stop = 0;
    int find_chromosome_flag = 1; // from chr1 or chr2 to find gene
    while(!stop){
        if (find_chromosome_flag){
            // when in chr1, find the gene that index w.r.t chr2
            next = chr2[find_gene(chr1,next)];

            chr2[find_gene(chr1,next)].Rbit = 0;

            find_chromosome_flag = 0;
        }else{
            // when in chr2, find the gene's index in chr1
            next = chr1[find_gene(chr1,next)];

            chr1[find_gene(chr1,next)].Rbit = 0;

            find_chromosome_flag = 1;
        }
        
        if (next.index == start_node.index){
            stop = 1;
        }
    }
}

int find_gene(trail chromosome, node gene){
// return the index where the gene at chromosome
    for (int i = 0;i < chromosome.size();i++){
        if (chromosome[i].index == gene.index){
            return i;
        } 
    }
    return -1;
}

void legalize_cx(trail &chr1,trail &chr2){
// if the Rbit == 1, exchange two gene between chr1 and chr2
    for (int i = 0;i < chr1.size();i++){
        if (chr1[i].Rbit == 1){
            node tmp = chr1[i];
            chr1[i] = chr2[i];
            chr2[i] = tmp;

            // set Rbit to 0;
            chr1[i].Rbit = 0;
            chr2[i].Rbit = 0;
        }
    }
}

void ox(trail &chr1, trail &chr2){
    // find the rest genes (Rbit = 1), and then perform legalize by gene's order 
    find_subtour(chr1,chr2);
    legalize_order(chr1,chr2);
}

void find_subtour(trail &chr1, trail &chr2){
    // init mark all Rbit to 1
    for (int i = 0;i < chr1.size();i++){
        chr1[i].Rbit = 1;
        chr2[i].Rbit = 1;
    }

    // randomly find the segment of the tour, set this subtour's Rbit = 0
    int effect_range = 7; // length of subtour
    int index_start = rand() % chr1.size();

    // make subtour Rbit = 0
    for (int i = 0;i < effect_range;i++){
        chr1[(index_start+i) % chr1.size()].Rbit = 0;
        chr2[(index_start+i) % chr2.size()].Rbit = 0;
    }
}

void legalize_order(trail &chr1, trail &chr2){
    // keep the information about chr1, chr2
    trail tmp_chr1 = chr1;
    trail tmp_chr2 = chr2;

    // perform the legalize to chr1 and chr2
    for (int i = 0;i < chr1.size();i++){
        if (chr1[i].Rbit == 1){
            // pick legal gene from chr2
            chr1[i] = LS_legalize_order(chr1,tmp_chr2);
            chr1[i].Rbit = 0;
        }

        if (chr2[i].Rbit == 1){
            // pick legal gene from chr1
            chr2[i] = LS_legalize_order(chr2,tmp_chr1);
            chr2[i].Rbit = 0;
        }
    }
}

node LS_legalize_order(trail illegal_chr, trail candidate_chr){
    bool isIllegal = 0;
    int index = 0;
    for (int i = 0;i < candidate_chr.size();i++){
        isIllegal = 0;
        for (int j = 0; j < illegal_chr.size();j++){
            if (illegal_chr[j].Rbit == 0 && candidate_chr[i].index == illegal_chr[j].index){
                isIllegal = 1;
            }
        }
        if(!isIllegal){
            index = i;
            return candidate_chr[i];
        }
    }
    return candidate_chr[index];
}

void mutate(trails &chromosomes){
    for (int i = 0;i < chromosomes.size();i++){
        // mutation rate = 0.2;
        if (rand() % 10 < 1){
            int const_i = rand() % chromosomes[i].size();
            int const_k = rand() % chromosomes[i].size();

            if (const_i == const_k){
                chromosomes[i] = TwoOptSwap(chromosomes[i],const_i,const_k+1);
            }else{
                chromosomes[i] = TwoOptSwap(chromosomes[i],const_i,const_k);
            }
        }
        TwoOpt(chromosomes[i]);
    }
}

void TwoOpt(trail &chromosome){
    int improve = 0;
    while (improve < 20){

        for (int start = 1;start < chromosome.size()-1;start++){
            for (int end = start + 1; end < chromosome.size()-1;end++){
                trail tmp_chr = TwoOptSwap(chromosome,start,end);

                if (evaluate_trail_distance(tmp_chr) < evaluate_trail_distance(chromosome)){
                    improve = 0;
                    chromosome = tmp_chr;
                }        
            }
        }
        improve++;
    }
}

trail TwoOptSwap(trail chromosome, const int& i, const int& k )
{

    int size = chromosome.size();
  
    trail new_chromosome = chromosome;

    // 2. take route[i] to route[k] and add them in reverse order to new_route
    int dec = 0;
    for ( int c = i; c <= k; c++ )
    {
        new_chromosome[c % new_chromosome.size()] = chromosome[ (k-dec) % chromosome.size()];
        dec++;
    }

    return new_chromosome;
  
}

trail pick_better_chromosome(trails pop_p){
// randomly return better individual in parent set
    double sum = 0.0;

    for (int i = 0;i < pop_p.size();i++){
        sum += evaluate_trail_distance(pop_p[i]);
    }

    double avg = sum / double(pop_p.size());

    

    std::vector<trail> candidate;
    for (int i = 0;i < pop_p.size();i++){
        if (evaluate_trail_distance(pop_p[i]) <= avg){
            candidate.push_back(pop_p[i]);
        }
    }

    return candidate[rand() % candidate.size()];
}

double evaluate_weight(node i, node j)
{
// computing the distance between two nodes
    return sqrt(pow(i.x - j.x, 2) + pow(i.y - j.y, 2));
}

double evaluate_trail_distance(trail t)
{
// computing the distance of the trail t
    t.push_back(t[0]);
    double distance = 0;
    for (int i = 0; i < t.size() - 1; i++)
    {
        distance += evaluate_weight(t[i], t[i + 1]);
    }
    return distance;
}

nodes readtsp(std::string path_to_file)
{
// read the dataset file
    std::ifstream input(path_to_file);
    std::string index;
    std::string x;
    std::string y;
    node n;

    nodes g;

    while (!input.eof())
    {
        input >> index;
        input >> x;
        input >> y;

        n.index = stod(index);
        n.x = stoi(x);
        n.y = stoi(y);
        n.Rbit = 0;

        g.push_back(n);
    }

    input.close();
    return g;
}

void outtsp(trail best)
{
    std::ofstream output("out.txt");
    for (int i = 0; i < best.size(); i++)
    {
        output << best[i].index << ' ';
        output << best[i].x << ' ';
        output << best[i].y << std::endl;
    }

    output << best[0].index << ' ';
    output << best[0].x << ' ';
    output << best[0].y << std::endl;

    output.close();
}

void showTrail(trail best)
{
// print the trail
    for (int i = 0; i < best.size(); i++)
    {
        std::cout << best[i].index << ' ';
    }
    std::cout << std::endl;
}

void showRbit(trail chromesome){
    for (int i = 0; i < chromesome.size(); i++)
    {
        std::cout << chromesome[i].Rbit << ' ';
    }
    std::cout << std::endl;

}

void plot(trail best)
{
// plot the trail by gnuplot
    outtsp(best);
    system("make plot");
}