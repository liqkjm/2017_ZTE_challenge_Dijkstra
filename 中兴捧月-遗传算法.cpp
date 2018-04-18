#include<stdio.h>
#include <stdlib.h> 
#include <memory.h>
 
#define MAXN 10000
#define INF 0x3f3f3f
#define SUM 10            //总共的染色体数量 
#define MAXloop 10000       //最大循环次数 
#define error 0.01        //若两次最优值之差小于此数则认为结果没有改变 
#define crossp 0.7        //交叉概率 
#define mp 0.1           //变异概率
 
int numnode = 0;
int dis[MAXN][MAXN];
int pathmatirx[MAXN][MAXN];
int numofnodepath[MAXN][MAXN];
int S ;   //起点
int T ;   //终点
int mustp[MAXN];
int nummustp = 0;   //必经点的个数，必经边算两个必经点
int NUMmustp = 0;   //个体中的染色体个数
 
int ret[1000];
int ptri = 0;
 
struct gen                        //定义染色体结构 
{ 
    int info[MAXN];               //染色体结构，表示从0开始经过这个染色体结构之后，到达T
    int cost;            //次染色体所对应的适应度函数值，在本题中为表达式的值 
    int numofnode;
};
struct gen gen_group[SUM];//定义一个含有20个染色体的组   
   
struct gen gen_result;     //记录次优的染色体
struct gen gen_result2;    //记录最优的染色体
int result_unchange_time; //记录在error前提下最优值为改变的循环次数 
 
struct mustedge{
    bool flag;
    int u;
    int v;
}ismustedge[MAXN];
 
 
//**************************************普通函数声明*****************************//
void reading();                            //读图
void floyd();                              //弗洛伊德求最短路径
void init();                               //初始化
int   randsign(float p);    //按照概率p产生随机数0、1，其值为1的概率为p  
int   randbit(int i,int j); //产生一个在i，j两个数之间的随机整数  
void findway( int a , int b );
int numofnode( int i , int j );
 
//**************************************遗传函数声明*****************************//
void gen_swap(gen *a,gen *b);                         //结构体交换
void gen_quicksort(gen *number,int left,int right);   //种群排序
void initiate();            //初始化函数，主要负责产生初始化种群  
void evaluation(int flag);  //评估种群中各染色体的适应度，并据此进行排序  
void Cross_group( gen &p, gen &q);               //交叉函数  
void selection();           //选择函数  
int  record();              //记录每次循环产生的最优解并判断是否终止循环  
void Varation_group(gen group[]);            //变异函数  
int Search_son( int path[], int len, int city);
int Search_son1( int path[], int len, int city);
void Rotate(int path[],int len, int m);
 
void evaluation(int flag) 
{ 
    int i , j , node1 , node2 ; 
    struct gen *genp;  
    genp = gen_group;
    for(i = 0 ; i < SUM ; i++)//计算各染色体对应的表达式值 
    { 
        genp[i].cost = 0;
        int cost = 0;
        int num = 0;
        for( j = 0 ; j < NUMmustp - 1 ; ++j ){
            if( !ismustedge[genp[i].info[j]].flag ){
                node1 = genp[i].info[j];
                node2 = genp[i].info[j + 1];
                cost += dis[node1][node2];
                num += numofnodepath[node1][node2];
            }
            else{
                node1 = genp[i].info[j];
                node2 = ismustedge[genp[i].info[j]].v;
                cost += dis[node1][node2];
                node1 = genp[i].info[j + 1];
                cost += dis[node2][node1];
                num += 1;
                num += numofnodepath[node2][node1];
            }
        }
        if( !ismustedge[genp[i].info[NUMmustp - 1]].flag ){
            node1 = genp[i].info[NUMmustp - 1];
            node2 = genp[i].info[0];
            cost += dis[node1][node2];
            num += numofnodepath[node1][node2];
        }
        else{
            node1 = genp[i].info[NUMmustp - 1];
            node2 = ismustedge[genp[i].info[NUMmustp - 1]].v;
            cost += dis[node1][node2];
            node1 = genp[i].info[0];
            cost += dis[node2][node1];
            num += 1;
            num += numofnodepath[node2][node1];
        }
        cost -= dis[T][S];
        genp[i].cost = cost;
        genp[i].numofnode = num;
    } 
    gen_quicksort(genp ,0 ,SUM-1 );              //对种群进行重新排序
} 
 
void calnodenum(){
    int i , j;
    for(i = 0 ; i < nummustp ; i++) {
        for(j = 0 ; j < nummustp; j++ ) {
            if(mustp[i] == mustp[j]){
                numofnodepath[mustp[i]][mustp[j]] = 0;
            }
            else{
                numofnodepath[mustp[i]][mustp[j]] = numofnode(mustp[i] , mustp[j]);
            }
            /*if(i == j){
                numofnodepath[i][j] = 0;
            }
            else{
                numofnodepath[i][j] = numofnode(i , j);
            }*/
        }
    }
}
 
int main(){
     
    int i , j;
    reading();
 
     
 
    floyd();
 
    calnodenum();
 
    result_unchange_time = 0;
    gen_result.cost = INF;
    gen_result2.cost = INF;
    initiate();
    evaluation( 0 );        //对初始化种群进行评估、排序  
    for( i = 0 ; i < MAXloop && result_unchange_time < 10 ; i++ ) 
    { 
        printf("第%d次迭代：",i);
        for(int  ii = 0 ; ii < SUM ; ++ii ){
            printf("染色体%d:    ",ii);
            for(j = 0 ; j < NUMmustp; ++j ){
                printf("%d ",gen_group[ii].info[j]);
            }
            printf("    花费为：%d，结点数为：%d\n",gen_group[ii].cost,gen_group[ii].numofnode);
        }printf("\n");
        float temp = 0;
        for (j = 0; j < SUM; j+= 1)  
        { 
            temp += gen_group[j].cost; 
        } 
 
        printf("本代平均花费为：%f\n",temp/SUM);
 
        printf("\n\n\n\n");
        if (gen_group[0].cost < gen_result.cost)  
        { 
            result_unchange_time = 0;
            memcpy(&gen_result, &gen_group[0], sizeof(gen)); 
        }
        else{
            result_unchange_time++;
        }
 
        for( j = 0; j < SUM ; ++j){
            if(gen_group[j].numofnode <= 9 && gen_group[j].cost < gen_result2.cost){
                result_unchange_time = 0;
                memcpy(&gen_result2, &gen_group[0], sizeof(gen));
            }
        }
 
        for (j = 0; j < SUM / 2; j+= 1)  
        { 
            Cross_group(gen_group[j], gen_group[ SUM - j -1]); 
        } 
        evaluation( 0 );
        Varation_group(gen_group);
        evaluation( 0 );
 
    }
 
    if(gen_result2.cost != INF){
        printf("有最优解：\n");
        memcpy(&gen_result, &gen_result2, sizeof(gen));
    }
    else{
        printf("无最优解，输出次优解：\n");
    }
 
    for(int ii=0;ii<NUMmustp - 1;ii++) 
    { 
        i = gen_result.info[ii];
        j = gen_result.info[ii+1];
        if(ismustedge[i].flag){
            //printf(",V%d,",i);
            ret[ptri++] = i;
            findway(ismustedge[i].v , j);
        }
        else
            findway(i , j);
    }
    i = gen_result.info[NUMmustp-1];
    j = gen_result.info[0];
    if(ismustedge[i].flag){
        //printf(",V%d,",i);
        ret[ptri++] = i;
        findway(ismustedge[i].v , j);
    }
    else
        findway(i , j);
    //printf("\n");printf("\n");
 
 
    int pos1 = Search_son1( ret, ptri, S);
    Rotate(ret, ptri, pos1);
    for( i = 0 ; i < ptri ; ++i){
        printf("%d->",ret[i]);
    }
    printf("\n");
    system("pause");
    return 0;
}
int numofnode( int i , int j ){
    int k , retnum = 0;
    k=pathmatirx[i][j];                               //取路径上Vi的后续Vk 
    if(k==-1) 
    { 
        printf("顶点%d 和 顶点%d 之间没有路径\n",i,j);//路径不存在  
    } 
    else 
    { 
        retnum++;
        while(k!=j) 
        {                         
            retnum++;
            k=pathmatirx[k][j];                   //求路径上下一顶点序号  
        } 
    }
    return retnum;
}
 
void findway( int i , int j ){
    int k ;
    k=pathmatirx[i][j];                               //取路径上Vi的后续Vk 
    if(k==-1) 
    { 
        printf("顶点%d 和 顶点%d 之间没有路径\n",i,j);//路径不存在  
    } 
    else 
    { 
        ret[ptri++] = i;
        while(k!=j) 
        {                         
            ret[ptri++] = k;
            k=pathmatirx[k][j];                   //求路径上下一顶点序号  
        } 
    } 
}
 
void Varation_group(gen group[]) 
{ 
    int i, j, k; 
    double temp; 
    //变异的数量，即，群体中的个体以PM的概率变异,变异概率不宜太大 
    int num = SUM * mp; 
   
    while (num--)  
    { 
        //确定发生变异的个体 
        k = rand() % SUM; 
   
        //确定发生变异的位 
        i = rand() % NUMmustp; 
        j = rand() % NUMmustp; 
   
        //exchange 
        temp  = group[k].info[i]; 
        group[k].info[i] = group[k].info[j];  
        group[k].info[j] = temp; 
    } 
} 
 
int Search_son1( int path[], int len, int city) 
{ 
    int i = 0; 
    for (i = 0; i < len; i++)  
    { 
        if (path[i] == city)  
        { 
            return i; 
        }
    } 
    return -1; 
} 
 
int Search_son( int path[], int len, int city) 
{ 
    int i = 0; 
    for (i = 0; i < len; i++)  
    { 
        if (path[i] == city)  
        { 
            return i; 
        }
        else if( ismustedge[ city ].flag && ismustedge[ city ].v == path[i] ){
            path[i] = ismustedge[ path[i] ].v;
            return i;
        }
    } 
    return -1; 
} 
 
//reverse a array 
//it's a auxiliary function for Rotate()  
void Reverse(int path[], int b, int e) 
{ 
    int temp; 
   
    while (b < e)  
    { 
        temp = path[b];  
        path[b] = path[e]; 
        path[e] = temp; 
   
        b++; 
        e--; 
    } 
} 
   
   
//旋转 m 位 
void Rotate(int path[],int len, int m) 
{ 
    if( m < 0 ) 
    { 
        return; 
    } 
    if (m > len)  
    { 
        m %= len; 
    } 
   
    Reverse(path, 0, m -1); 
    Reverse(path, m, len -1); 
    Reverse(path, 0, len -1); 
} 
 
void Cross_group( gen &p, gen &q) 
{ 
    //for(int ii = 0 ; ii < NUMmustp ; ++ ii){
    //  printf("%d ",p.info[ii]);
    //}printf("\n");
    //for(int ii = 0 ; ii < NUMmustp ; ++ ii){
    //  printf("%d ",q.info[ii]);
    //}printf("\n");printf("\n");printf("\n");
 
    int i = 0; 
    int pos1, pos2; 
    int len = NUMmustp; 
    int first; 
 
    double len1 ,len2 ;
    if( ismustedge[ p.info[0] ].flag )
        len1 = dis[ismustedge[ p.info[0]].v][ p.info[1] ]; 
    else
        len1 = dis[p.info[0] ][ p.info[1] ]; 
    if( ismustedge[ q.info[0] ].flag )
        len2 = dis[ismustedge[ q.info[0]].v][ q.info[1] ];
    else
        len2 = dis[q.info[0] ][ q.info[1] ]; 
 
    if (len1 <= len2)  
    { 
        first = p.info[0]; 
    } 
    else 
    { 
        first = q.info[0]; 
    } 
    pos1 = Search_son( p.info + i, len, first); 
    pos2 = Search_son( q.info + i, len, first); 
   
    Rotate(p.info + i, len, pos1); 
    Rotate(q.info + i, len, pos2); 
 
    while ( --len > 1)  
    {
        i++; 
        int span1 , span2 ;
 
        int temp;
        if(ismustedge[ p.info[i - 1] ].flag){
            temp = ismustedge[ p.info[i - 1] ].v;
        }
        else{
            temp = p.info[i - 1];
        }
        span1 = dis[temp][ p.info[i] ]; 
 
        if(ismustedge[ q.info[i - 1] ].flag){
            temp = ismustedge[ q.info[i - 1] ].v;
        }
        else{
            temp = q.info[i - 1];
        }
        span2 = dis[temp][ q.info[i] ]; 
 
        if ( span1 <= span2 ) 
        { 
            pos2 = Search_son( q.info + i, len, p.info[i]); 
            Rotate(q.info + i, len, pos2); 
        } 
        else 
        { 
            pos1 = Search_son( p.info + i, len, q.info[i]); 
            Rotate(p.info + i, len, pos1); 
        }
    } 
     
    Rotate(q.info, NUMmustp, rand() % NUMmustp);
} 
  
void initiate(){
    /*********第一个初始解定义为按离起点的距离的顺序*********************/
    bool *flag = NULL;                                 //用于标记必经点中的某点是否被选进染色体中
    flag = (bool*)malloc(sizeof(bool) * nummustp);
    if(flag == NULL){
        printf("error initiate\n"); 
        exit(1); 
    }
    for(int i = 0 ; i < nummustp; ++i )flag[i] = false;
    int i , j , z ;
    gen_group[0].info[NUMmustp - 1] = T ; flag[nummustp - 1] = true; flag[0] = true;
    for( i = 0 ; i < NUMmustp - 1 ; ++i ){
        int min = INF;int k = INF;
        for( j = 0 ; j < nummustp ; ++j ){
            if(!flag[j] && dis[S][mustp[j]] < min ){
                min = dis[S][mustp[j]];
                k = j;
            }
        }
         
        if(k != INF){
            if( !ismustedge[mustp[k]].flag ){     //如果是必经点
                gen_group[0].info[i] = mustp[k];
                flag[k] = true;
            }
            else{                                 //如果是必经边
                gen_group[0].info[i] = mustp[k];
                flag[k] = true;
                for(z = 0 ; z < nummustp ; ++z ){
                    if( mustp[z] == ismustedge[mustp[k]].v ){
                        flag[z] = true;break;
                    }
                }
            }
        }
        else{
            break;
        }
    }
 
    /*********随机生成剩余初始解*********************/
    int k;
    for( i = 1 ; i < SUM ; ++i ){
        for(j = 0 ; j < nummustp; ++j )flag[j] = false;
        gen_group[i].info[NUMmustp - 1] = T; flag[0] = true; flag[nummustp - 1] = true;
        for(j = 0 ; j < NUMmustp - 1 ; ++j ){
            k = randbit(0 , nummustp-1);
            while( flag[k] ){
                k = randbit(0 , nummustp-1);
            }
 
            if( !ismustedge[mustp[k]].flag ){     //如果是必经点
                gen_group[i].info[j] = mustp[k];
                flag[k] = true;
            }
            else{                                 //如果是必经边
                gen_group[i].info[j] = mustp[k];
                flag[k] = true;
                for(z = 0 ; z < nummustp ; ++z ){
                    if( mustp[z] == ismustedge[mustp[k]].v ){
                        flag[z] = true;break;
                    }
                }
            }
        }
    }
 
    free(flag);
    flag = NULL;
}
 
void reading(){
    FILE *fp; 
    if(NULL == (fp = fopen("case0.txt", "r"))) 
    { 
        printf("error reading\n"); 
        exit(1); 
    } 
 
    S = 0; numnode = 0;
    char ch; 
    while( '\n' != (ch=fgetc(fp)) )  //总的点的个数
    { 
        numnode = numnode*10 + ch - '0'; 
        //printf("%c", ch); 
    } 
    T = numnode - 1;
    //printf("起点为%d，终点为%d，结点数目为%d\n" , S, T, numnode);
   
    ch=fgetc(fp);
 
    nummustp = 0;
    memset(mustp,0,sizeof(mustp));
    mustp[++nummustp] = S;++NUMmustp;
    while( '\n' != (ch=fgetc(fp)) )    //读取必经点
    { 
        if(ch == ' '){
            ++nummustp;++NUMmustp;
        }
        else{
            mustp[nummustp] = mustp[nummustp]*10 + ch - '0'; 
            //printf("%c", ch); 
        }
    }
 
    ch=fgetc(fp);
 
    init();
    int temp[3] = {0,0,0} , j = 0;
    while( '\n' != (ch=fgetc(fp)) )    //读取图
    { 
        temp[0] = 0 ,temp[1] = 0 ,temp[2] = 0 , j = 0;
        while(ch != '\n'){
            if( ch == ' ' ){
                ++j;
            }
            else{
                temp[j] = temp[j]*10 + ch - '0';
            }
            ch = fgetc(fp);
        }
        dis[temp[0]][temp[1]] = temp[2];
        dis[temp[1]][temp[0]] = temp[2];
        pathmatirx[temp[0]][temp[1]] = temp[0];
        pathmatirx[temp[1]][temp[0]] = temp[1];
    } 
 
    while( '\n' != (ch=fgetc(fp)) )   //必经边的权值设为0
    { 
        temp[0] = 0 ,temp[1] = 0 ,temp[2] = 0 , j = 0;
        while(ch != '\n'){
            if( ch == ' ' ){
                ++j;
            }
            else{
                temp[j] = temp[j]*10 + ch - '0';
            }
            ch = fgetc(fp);
        }
        mustp[++nummustp] = temp[0];
        mustp[++nummustp] = temp[1];
        ++NUMmustp;
        ismustedge[temp[0]].flag = true;
        ismustedge[temp[0]].u = temp[0];
        ismustedge[temp[0]].v = temp[1];
        ismustedge[temp[1]].flag = true;
        ismustedge[temp[1]].u = temp[1];
        ismustedge[temp[1]].v = temp[0];
    } 
 
    while( '\n' != (ch=fgetc(fp)) )   //惩罚边的权值设为INF
    { 
        temp[0] = 0 ,temp[1] = 0 ,temp[2] = 0 , j = 0;
        while(ch != '\n'){
            if( ch == ' ' ){
                ++j;
            }
            else{
                temp[j] = temp[j]*10 + ch - '0';
            }
            ch = fgetc(fp);
        }
        dis[temp[0]][temp[1]] = INF;
        dis[temp[1]][temp[0]] = INF;
    }
    ++NUMmustp;
    mustp[++nummustp] = T;
    ismustedge[S].flag = true;
    ismustedge[S].u = S;
    ismustedge[S].v = T;
    ismustedge[T].flag = true;
    ismustedge[T].u = T;
    ismustedge[T].v = S;
    ++nummustp;            //必经点的数量+1
    //printf("\n");
    fclose(fp);
}
 
void floyd(){                                  //计算所有顶点到所有顶点的最短路径
    int k , i , j; 
 
    for( i = 0 ; i <= numnode ; ++i ){
        for( j = 0 ; j <= numnode ; ++j ){
            if(dis[i][j] < INF)
                pathmatirx[i][j] = j;               //初始化路径数组
            else
                pathmatirx[i][j] = -1;
        }
    }
 
    for(k = 0 ; k < numnode ; ++k ) 
    { 
        for( i = 0 ; i < numnode ; ++i ) 
        { 
            for( j = 0 ; j < numnode ; ++j ) 
            { 
                if( dis[i][j] > dis[i][k] + dis[k][j] ) { 
                    //如果经过下标为k顶点路径比原两点间路径更短
                    dis[i][j] = dis[i][k] + dis[k][j];   //更新当前两点间权值
                    pathmatirx[i][j] = pathmatirx[i][k]; //路径设置为经过下标为K的顶点
                }
            } 
        } 
    } 
} 
 
void gen_swap(gen *a,gen *b)
{
    gen temp;
    temp = *a;
    *a = *b;
    *b = temp;
}
 
void gen_quicksort(gen *number,int left,int right)//快速排序，用于结构体排序
{
    int i ,j ,s ;
    if(left < right)
    {
        s = number[(left + right) / 2].cost ,j = right+1 ,i = left-1;
        while(1)
        {
            while(number[++i].cost < s )  ;
            while(number[--j].cost > s )  ;
            if(i>=j)
                break;
            gen_swap(&number[i] ,&number[j] );
        }
        gen_quicksort(number ,left ,i-1 );
        gen_quicksort(number ,j+1 ,right );
    }
}
 
void init(){
    int i , j ;
    for( i = 0 ; i <= numnode ; ++i ){
        for( j = 0 ; j <= numnode ; ++j ){
            //pathmatirx[i][j] = j;               //初始化路径数组
            if(i == j ) dis[i][j] = 0;
            else dis[i][j] = INF;
        }
    }
    for(i = 0 ; i < MAXN ; ++i ){
        ismustedge[i].flag = false;
    }
}
 
int randsign(float p)//按概率p返回1 
{ 
    if(rand() > (p * 32768)) 
        return 0; 
    else return 1; 
} 
int randbit(int i, int j)//产生在i与j之间的一个随机数 
{ 
    int a , l; 
    l = j - i + 1; 
    a = i + rand() * l / 32768; 
    return a; 
}
