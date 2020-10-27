#include "classicrbn.h"

void initModular(RBN* net, int nm, float k, float l, int m)
{
    register int i;
    net->nm         = nm;
    net->k          = k;
    net->l          = l;
    net->m          = m;
    net->n          = nm * m;
    net->p          = 0.5f;
    net->conn       = ( char** ) calloc( net->n, sizeof( char* ) );
    net->rules      = ( char** ) calloc( net->n, sizeof( char* ) );
    net->state      = ( char*  ) calloc( net->n, sizeof( char  ) );
    net->nextState  = ( char*  ) calloc( net->n, sizeof( char  ) );
    net->links      = ( int**  ) calloc( net->n, sizeof( int*  ) );
    for (i = 0; i < net->n; i++)
    {
        net->conn[i] = ( char* ) calloc( net->n, sizeof( char  ) );
        net->links[i]= ( int* ) calloc( net->k + 1, sizeof( int ) );//TODO: check type of links (int) or char 9abr2013
    }
    initM(net, 0, net->n);
}

void init(RBN* net, int n, float k)
{
    initModular(net, n, k, 0.0f, 1);
}

void initM(RBN* net, const int from, const int to)
{
    register int i, j;
    for (i = from; i < to; i++)
    {
        net->state[i] = net->nextState[i] = 0;
        for (j = 0; j < to; j++)
        {
            net->conn[i][j] = 0;
            net->conn[j][i] = 0;
        }
    }
}

void initRules(RBN* net, float prob)
{
    register int i, j;
    for (i = 0; i < net->n; i++)
    {
        net->rules[i] = ( char* ) realloc(net->rules[i], pow(2, kAt(net, i)) * sizeof( char ) );
        for (j = 0; j < pow(2, kAt(net, i)); j++)
        {
            net->rules[i][j] = ((float)rand() / RAND_MAX) < prob;
        }
    }
}

void randomInitial(RBN* net)
{
    register int i;
    for (i = 0; i < net->n; i++)
    {
        net->state[i] = rand() % 2;
    }
}

void step(RBN* net)
{
    register int i, j, ruleno, nodeno;
    for (i = 0; i < net->n; i++)
    {
        ruleno = 0;
        nodeno = 1;
        for (j = net->n - 1; j > -1; j--)
        {
            if ( net->conn[i][j] )
            {
                if (net->state[j] == 1)
                {
                    ruleno += nodeno;
                }
                //ruleno += net->state[j] * nodeno;
                nodeno *= 2;
            }
        }
        net->nextState[i] = net->rules[i][ruleno];//apply rule
    }
    memcpy(net->state, net->nextState, sizeof( char ) * net->n );
}

/**
 * Adyacent list
 */
void next(RBN* net)
{
    register int i, j, ruleno, nodeno;
    for (i = 0; i < net->n; i++)
    {
        ruleno = 0;
        nodeno = 1;
        for (j = 0; net->links[i][j] != -1; j++)
        {
            ruleno += net->state[net->links[i][j]] * nodeno;
            nodeno *= 2;
        }
        net->nextState[i] = net->rules[i][ruleno];//apply rule
    }
    memcpy(net->state, net->nextState, sizeof( char ) * net->n );
}

void generateLinksAdyacentList(RBN* net)
{
    register int i, j, k;
    for (i = 0; i < net->n; i++)
    {
        net->links[i] = ( int* ) realloc(net->links[i], (kAt(net, i) + 1) * sizeof( int ) );
        for (j = net->n - 1, k = 0; j > -1; j--)
        {
            if ( net->conn[i][j] )
            {
                net->links[i][k++] = j;
            }
        }
        net->links[i][k] = -1;
    }
    /*register int i, j;
    for (i = 0; i < net->n; i++)
    {
        net->links[i] = ( char** ) calloc(kAt(net, i) + 1, sizeof( char* ) );
        for (j = 0; j < net->n; j++)
        {
            if ( net->conn[i][j] )
            {
                net->links[i][j] = net->state + j;
            }
        }
        net->links[i][j + 1] = NULL;
    }*/
}

/**
 * Randomly generate connections for a {@link redrbn.ClassicRBN.Topology#HOMOGENEOUS}
 * topology network.
 * <p>Connects exactly k nodes to each node.</p>
 */
void generateHomogeneousConnectionMatrix(RBN* net)
{
    register int i, j, tmp;
    for (i = 0; i < net->n; i++)
    {
        j = 0;
        while (j < net->k)
        {
            tmp = rand() % net->n; //get random node
            if (!net->conn[i][tmp]) //let's see if we haven't used it
            {
                net->conn[i][tmp] = 1;
                j++;
            }
        }
    }
}

/**
 * Randomly generates connections for a {@link redrbn.ClassicRBN.Topology#NORMAL} topology network,
 * connecting <b>on average</b> <i>k</i> nodes to each node.
 */
void generateNormalConnectionMatrix(RBN* net)
{
    register int i;
    register int j;
    for (i = 0; i < net->n; i++)
    {
        for (j = 0; j < net->n; j++)
        {
            net->conn[i][j] = ( ((float)rand() / RAND_MAX) * net->n < net->k);
        }
    }
}

/**
 * Randomly generate connections for a {@link redrbn.ClassicRBN.Topology#SCALE_FREE}
 * topology network.
 * <p>Generates scale-free topology for inputs, <i>k</i> stands for <i>&#x3B3;</i>.</p>s
 */
void generateScaleFreeConnectionMatrix(RBN* net)
{
    register int i, j, tmp;
    int Ks[net->n];//number of connections that each node will receive
    for (i = 0; i < net->n; i++)
    {
        Ks[i] = (int) floor(pow( ((float)rand() / RAND_MAX), -(1 / (net->k - 1))));
        if (Ks[i] > net->n)
        {
            Ks[i] = net->n; //bound number of connections to N (only one K per node allowed...)
        }
    }
    for (i = 0; i < net->n; i++)
    {
        j = 0;
        while (j < floor(Ks[i]))
        {
            tmp = (int) floor( ((float)rand() / RAND_MAX) * net->n); //get random node
            if (!net->conn[i][tmp]) //let's see if we haven't used it
            {
                net->conn[i][tmp] = 1;
                j++;
            }
        }
    }
}

/**
 * Randomly generates the connection matrix for a network whose topology is coupled.
 */
void generateCoupledConnectionMatrix(RBN* net)
{
    int i, j, x, y;
    if (net->m == 1)
    {
        generateNormalConnectionMatrix(net);
    }
    else
    {
        int complemento = (net->m - 1) * net->nm;// todos aquellos que no son de Nets sino de interconexiones.
        /*almacenamos aquí los elejibles, que son repetibles por cada módulo m_i, n veces, dado que es información redundante tener un arreglo de M*n x complemento*/
        int** libres = (int**) calloc( net->n, sizeof(int*) );
        for (i = 0; i < net->n; i++)
        {
            libres[i] = (int*) calloc( complemento, sizeof(int) );
        }
        initMatrix(libres, net->n, complemento, -1);
        int libreContador = 0;
        for (x = 0, y = 0; x < net->m; x++) //modules
        {
            for (i = 0; i < net->n; i++) //rows
            {
                for (j = 0; j < net->n; j++) //columns
                {
                    if (x == y) //it should be an inner-matrix (i.e. a matrix in the diagonal) and
                    {
                        //a number between 0 and n
                        net->conn[i][j] = ( ((float)rand() / RAND_MAX) * net->nm < net->k);
                    }
                    else
                    {
                        /* Generate intermodule connectivity matrix to later on choose L connections*/
                        libres[i][libreContador] = j;//mark it as free
                        libreContador++;
                        if (libreContador % complemento == 0)
                        {
                            libreContador = 0;
                        }
                    }
                    if ((j + 1) % net->nm == 0)
                    {
                        y++;
                    }
                }
                if ((i + 1) % net->nm == 0)
                {
                    x++;
                }
                y = 0;
            }
        }
        generateInterModuleConnections(net, libres);
    }
}

/**
 * Connects the each module with its corresponding L module-neighbours.
 *
 * @param freePositions an array containing the free positions of the
 *          connectivity matrix, that is the positions in the matrix where
 *          a value of <i>1</i> would represent a link with another module.
 */
void generateInterModuleConnections(RBN* net, int** freePositions)
{
    int complement = (net->m - 1) * net->nm;
    int modCont = 0;
    int x;
    register int i, j;
    for (i = 0; i < net->n; i++)
    {
        for (j = 0; j < complement; j++)
        {
            do
            {
                x = rand() * net->nm + modCont;// TODO check rand(nm)
            }
            while ( freePositions[x][j] == -1 );
            net->conn[x][freePositions[x][j]] = ( ((float)rand() / RAND_MAX) * net->nm * complement < net->l);
            freePositions[x][j] = -1;//just tu mark it, though it is not really necessary
        }
        if ((i + 1) % net->nm == 0)
        {
            modCont += net->nm;
        }
    }
}

/**
 * Gets an array with the data indexes shuffled of the input samples
 * @return cards indexes randombly ordered
 */
void shuffle(int* cards, int n)
{
    register int i, r, temp;
    for (i = 0; i < n; i++)
    {
        r = rand() % n;
        temp = cards[i];
        cards[i] = cards[r];
        cards[r] = temp;
    }
}

/**
 * Duplicate functionality of a node adding the redundant node at least.
 * @param net The network that will be added a node.
 * @param deg Degree of Degeneracy for \link {Redundancy}
 */
void addNode(RBN* net, float* deg)
{
    register int i, j, k, aux, size, degNode = rand() % net->n;
    int* ones   = ( int* ) calloc( net->n, sizeof( int ) );
    int* zero   = ( int* ) calloc( net->n, sizeof( int ) );
    for (i = 0; i < net->n; i++)
    {
        net->conn[net->n][i] = net->conn[degNode][i];
        net->conn[i][net->n] = net->conn[i][degNode];
    }
    net->conn[net->n][net->n] = net->conn[degNode][degNode];

    for (i = j = k = 0; i < net->n; i++)
    {
        if ( net->conn[net->n][i] )
        {
            ones[j++] = i;
        }
        else
        {
            zero[k++] = i;
        }
    }
    shuffle(ones, j);
    shuffle(zero, k);
    aux = (j < k ? j : k) * (1 - deg[INPUT_DEGENERACY]);
    for (i = 0; i < aux; i++)
    {
        net->conn[net->n][ones[i]] = 0;
        net->conn[net->n][zero[i]] = 1;
    }

    for (i = 0; i < net->n; i++)
    {
        ones[i] = 0;
        zero[i] = 0;
    }

    for (i = j = k = 0; i < net->n; i++)
    {
        if (net->conn[i][net->n])
        {
            ones[j++] = i;
        }
        else
        {
            zero[k++] = i;
        }
    }
    shuffle(ones, j);
    shuffle(zero, k);
    aux = (j < k ? j : k) * (1 - deg[OUTPUT_DEGENERACY]);
    for (i = 0; i < aux; i++)
    {
        net->conn[ones[i]][net->n] = 0;
        net->conn[zero[i]][net->n] = 1;
    }

    for (i = 0; i < net->n; i++) // Copy function deg
    {
        if ( net->conn[i][net->n] )
        {
            aux = pow(2, kAt(net, i));
            net->rules[i] = (char*)realloc(net->rules[i], (aux * 2) * sizeof(char));
            for (j = 0; j < aux; j++)
            {
                net->rules[i][j + aux] = net->rules[i][j];
            }
        }
    }
    if (net->conn[net->n][net->n])
    {
        size = pow(2, kAt(net, net->n) + 1);
    }
    else
    {
        size = pow(2, kAt(net, net->n));
    }
    net->rules[net->n] = (char*)realloc(net->rules[net->n], size * sizeof(char));
    aux = pow(2, kAt(net, degNode));
    for (i = 0; i < size; i++)
    {
        net->rules[net->n][i] = net->rules[degNode][i % aux];
    }

    ones = (int*) realloc( ones, pow(2, kAt(net, net->n)) * sizeof(int) );
    for (i = j = k = 0; i < pow(2, kAt(net, net->n)); i++)
    {
        ones[j++] = i;
    }
    shuffle(ones, j);
    aux = j * (1 - deg[FUNCTION_DEGENERACY]);
    for (i = 0; i < aux; i++)
    {
        net->rules[net->n][ones[i]] = !net->rules[net->n][ones[i]];
    }

    free(ones);
    free(zero);
    net->n = net->nm = net->n + 1;
}

void toString(RBN* net)
{
    register int i, j;
    printf("  N: %d, K: %.2f, M: %d, L: %.2f, P: %.2f, Q:.2f\n\n", net->n, net->k, net->m, net->l, net->p);
    printf("state: ");
    for (i = 0; i < net->n; i++)
    {
        printf("%d",    net->state[i]);
    }
    printf("\nnet:\n");
    for (i = 0; i < net->n; i++)
    {
        for (j = 0; j < net->n; j++)
        {
            printf("%d", net->conn[i][j]);
        }
        printf("\n");
    }
    printf("rules:\n");
    for (i = 0; i < net->n; i++)
    {
        for (j = 0; j < pow(2, kAt(net, i)); j++)
        {
            printf("%d", net->rules[i][j]);
        }
        printf("\n");
    }
}

int getK(RBN* net)
{
    register int i, j, cont = 0;
    for (i = 0; i < net->n; i++)
    {
        for (j = 0; j < net->n; j++)
        {
            cont += net->conn[i][j];
        }
    }
    return cont;
}

float kTot(RBN* net)
{
    return net->k + net->l / net->n;
}

int kAt(RBN* net, int node)
{
    register int i;
    int count = 0;
    for (i = 0; i < net->n; i++)
    {
        if (net->conn[node][i])
        {
            count++;
        }
    }
    return count;
}

/**
 * @return true if node has connections
 */
_Bool anyKAt(RBN* net, int node)
{
    register int j;
    for (j = 0; j < net->n; j++)
    {
        if (net->conn[node][j])
        {
            return 1;
        }
    }
    return 0;
}

/**
 * @return column of nth node
 */
int kTh(RBN* net, int row, int nth)
{
    register int j, count = 0;
    for (j = 0; j < net->n; j++)
    {
        if (net->conn[row][j])
        {
            count++;
        }
        if (nth < count)
        {
            return j;
        }
    }
    return 0;
}

/**
 * Fills the array which contains the "free" spaces on the connection matrix with the value <i>-1</i>.
 *
 * @param matrix The matrix to fill.
 */
void initMatrix(int** matrix, int rows, int cols, int value)
{
    register int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            matrix[i][j] = value;
        }
    }
}

void print(char* value, int n )
{
    register int i;
    for (i = 0; i < n; i++)
    {
        printf("%d", value[i]);
    }
    printf("\n");
}

void printMxN(char** value, int rows, int cols )
{
    register int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            printf("%d", value[i][j]);
        }
        printf("\n");
    }
}

/**
 * Generate Redundant Matrix
 */
char** calcRed(RBN* net)
{
    register int i, i2, parts, size;
    int inputK;//keeps number of inputs
    int redundant;//flag to check if link is redundant
    int tmp1, tmp2;
    char** red = (char**) calloc(net->n, sizeof( char* ) );
    for (i = 0; i < net->n; i++)  //init matrix
    {
        red[i] = (char*) calloc(net->n, sizeof( char ) );
    }
    for (i = 0; i < net->n; i++)
    {
        inputK = kAt(net, i);
        if (anyKAt(net, i))  //if any link, check if they are redundant
        {
            for (i2 = 0; i2 < inputK; i2++)  //browse all inputs
            {
                redundant = 1;
                for (parts = 0; parts < (int) pow(2, i2); parts++)
                {
                    for (size = 0; size < (int) pow(2, inputK - i2 - 1); size++)
                    {
                        tmp1 = net->rules[i][size + (int) (pow(2, inputK - i2) * parts)];
                        tmp2 = net->rules[i][size + (int) (pow(2, inputK - i2) * parts)
                                             + (int) pow(2, inputK - i2 - 1)];
                        if (tmp1 != tmp2)
                        {
                            redundant = 0;
                        }
                    }
                }
                if (redundant)  //mark in the redundancy matrix
                {
                    red[i][kTh(net, i, i2)] = 1;
                }
            }
        }
    }
    return red;
}

void removeRedLink(RBN* net, char** red)
{
    /*register int i, j, parts, size;
    int nodeTBR;//node to be removed
    if (countRed() > 0)
    {
        nodeTBR = rand() % countRed(net);
        int count = 0;
        for (i = 0; i < net->n; i++)
        {
            for (j = 0; j < net->n; j++)
            {
                if (red[i][j])
                {
                    if (count == nodeTBR)  //remove in matrix and rules...
                    {
                        int inputK = kAt(net, i);
                        int inputNo = inputNo(i, j);
                        for (parts = (int) pow(2, inputNo) - 1; parts > -1; parts--)
                        {
                            for (size = (int) pow(2, inputK - inputNo - 1) - 1; size > -1; size--)
                            {
                                net->rules[i].remove(size + (int) (pow(2, inputK - inputNo) * parts));
                            }
                        }
                        net->conn[i][j] = 0;
                    }
                    count++;
                }
            }
        }
    }*/
}
