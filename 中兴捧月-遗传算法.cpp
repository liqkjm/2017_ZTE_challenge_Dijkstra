#include<stdio.h>
#include <stdlib.h> 
#include <memory.h>
 
#define MAXN 10000
#define INF 0x3f3f3f
#define SUM 10            //�ܹ���Ⱦɫ������ 
#define MAXloop 10000       //���ѭ������ 
#define error 0.01        //����������ֵ֮��С�ڴ�������Ϊ���û�иı� 
#define crossp 0.7        //������� 
#define mp 0.1           //�������
 
int numnode = 0;
int dis[MAXN][MAXN];
int pathmatirx[MAXN][MAXN];
int numofnodepath[MAXN][MAXN];
int S ;   //���
int T ;   //�յ�
int mustp[MAXN];
int nummustp = 0;   //�ؾ���ĸ������ؾ����������ؾ���
int NUMmustp = 0;   //�����е�Ⱦɫ�����
 
int ret[1000];
int ptri = 0;
 
struct gen                        //����Ⱦɫ��ṹ 
{ 
    int info[MAXN];               //Ⱦɫ��ṹ����ʾ��0��ʼ�������Ⱦɫ��ṹ֮�󣬵���T
    int cost;            //��Ⱦɫ������Ӧ����Ӧ�Ⱥ���ֵ���ڱ�����Ϊ���ʽ��ֵ 
    int numofnode;
};
struct gen gen_group[SUM];//����һ������20��Ⱦɫ�����   
   
struct gen gen_result;     //��¼���ŵ�Ⱦɫ��
struct gen gen_result2;    //��¼���ŵ�Ⱦɫ��
int result_unchange_time; //��¼��errorǰ��������ֵΪ�ı��ѭ������ 
 
struct mustedge{
    bool flag;
    int u;
    int v;
}ismustedge[MAXN];
 
 
//**************************************��ͨ��������*****************************//
void reading();                            //��ͼ
void floyd();                              //�������������·��
void init();                               //��ʼ��
int   randsign(float p);    //���ո���p���������0��1����ֵΪ1�ĸ���Ϊp  
int   randbit(int i,int j); //����һ����i��j������֮����������  
void findway( int a , int b );
int numofnode( int i , int j );
 
//**************************************�Ŵ���������*****************************//
void gen_swap(gen *a,gen *b);                         //�ṹ�彻��
void gen_quicksort(gen *number,int left,int right);   //��Ⱥ����
void initiate();            //��ʼ����������Ҫ���������ʼ����Ⱥ  
void evaluation(int flag);  //������Ⱥ�и�Ⱦɫ�����Ӧ�ȣ����ݴ˽�������  
void Cross_group( gen &p, gen &q);               //���溯��  
void selection();           //ѡ����  
int  record();              //��¼ÿ��ѭ�����������ŽⲢ�ж��Ƿ���ֹѭ��  
void Varation_group(gen group[]);            //���캯��  
int Search_son( int path[], int len, int city);
int Search_son1( int path[], int len, int city);
void Rotate(int path[],int len, int m);
 
void evaluation(int flag) 
{ 
    int i , j , node1 , node2 ; 
    struct gen *genp;  
    genp = gen_group;
    for(i = 0 ; i < SUM ; i++)//�����Ⱦɫ���Ӧ�ı��ʽֵ 
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
    gen_quicksort(genp ,0 ,SUM-1 );              //����Ⱥ������������
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
    evaluation( 0 );        //�Գ�ʼ����Ⱥ��������������  
    for( i = 0 ; i < MAXloop && result_unchange_time < 10 ; i++ ) 
    { 
        printf("��%d�ε�����",i);
        for(int  ii = 0 ; ii < SUM ; ++ii ){
            printf("Ⱦɫ��%d:    ",ii);
            for(j = 0 ; j < NUMmustp; ++j ){
                printf("%d ",gen_group[ii].info[j]);
            }
            printf("    ����Ϊ��%d�������Ϊ��%d\n",gen_group[ii].cost,gen_group[ii].numofnode);
        }printf("\n");
        float temp = 0;
        for (j = 0; j < SUM; j+= 1)  
        { 
            temp += gen_group[j].cost; 
        } 
 
        printf("����ƽ������Ϊ��%f\n",temp/SUM);
 
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
        printf("�����Ž⣺\n");
        memcpy(&gen_result, &gen_result2, sizeof(gen));
    }
    else{
        printf("�����Ž⣬������Ž⣺\n");
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
    k=pathmatirx[i][j];                               //ȡ·����Vi�ĺ���Vk 
    if(k==-1) 
    { 
        printf("����%d �� ����%d ֮��û��·��\n",i,j);//·��������  
    } 
    else 
    { 
        retnum++;
        while(k!=j) 
        {                         
            retnum++;
            k=pathmatirx[k][j];                   //��·������һ�������  
        } 
    }
    return retnum;
}
 
void findway( int i , int j ){
    int k ;
    k=pathmatirx[i][j];                               //ȡ·����Vi�ĺ���Vk 
    if(k==-1) 
    { 
        printf("����%d �� ����%d ֮��û��·��\n",i,j);//·��������  
    } 
    else 
    { 
        ret[ptri++] = i;
        while(k!=j) 
        {                         
            ret[ptri++] = k;
            k=pathmatirx[k][j];                   //��·������һ�������  
        } 
    } 
}
 
void Varation_group(gen group[]) 
{ 
    int i, j, k; 
    double temp; 
    //���������������Ⱥ���еĸ�����PM�ĸ��ʱ���,������ʲ���̫�� 
    int num = SUM * mp; 
   
    while (num--)  
    { 
        //ȷ����������ĸ��� 
        k = rand() % SUM; 
   
        //ȷ�����������λ 
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
   
   
//��ת m λ 
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
    /*********��һ����ʼ�ⶨ��Ϊ�������ľ����˳��*********************/
    bool *flag = NULL;                                 //���ڱ�Ǳؾ����е�ĳ���Ƿ�ѡ��Ⱦɫ����
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
            if( !ismustedge[mustp[k]].flag ){     //����Ǳؾ���
                gen_group[0].info[i] = mustp[k];
                flag[k] = true;
            }
            else{                                 //����Ǳؾ���
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
 
    /*********�������ʣ���ʼ��*********************/
    int k;
    for( i = 1 ; i < SUM ; ++i ){
        for(j = 0 ; j < nummustp; ++j )flag[j] = false;
        gen_group[i].info[NUMmustp - 1] = T; flag[0] = true; flag[nummustp - 1] = true;
        for(j = 0 ; j < NUMmustp - 1 ; ++j ){
            k = randbit(0 , nummustp-1);
            while( flag[k] ){
                k = randbit(0 , nummustp-1);
            }
 
            if( !ismustedge[mustp[k]].flag ){     //����Ǳؾ���
                gen_group[i].info[j] = mustp[k];
                flag[k] = true;
            }
            else{                                 //����Ǳؾ���
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
    while( '\n' != (ch=fgetc(fp)) )  //�ܵĵ�ĸ���
    { 
        numnode = numnode*10 + ch - '0'; 
        //printf("%c", ch); 
    } 
    T = numnode - 1;
    //printf("���Ϊ%d���յ�Ϊ%d�������ĿΪ%d\n" , S, T, numnode);
   
    ch=fgetc(fp);
 
    nummustp = 0;
    memset(mustp,0,sizeof(mustp));
    mustp[++nummustp] = S;++NUMmustp;
    while( '\n' != (ch=fgetc(fp)) )    //��ȡ�ؾ���
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
    while( '\n' != (ch=fgetc(fp)) )    //��ȡͼ
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
 
    while( '\n' != (ch=fgetc(fp)) )   //�ؾ��ߵ�Ȩֵ��Ϊ0
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
 
    while( '\n' != (ch=fgetc(fp)) )   //�ͷ��ߵ�Ȩֵ��ΪINF
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
    ++nummustp;            //�ؾ��������+1
    //printf("\n");
    fclose(fp);
}
 
void floyd(){                                  //�������ж��㵽���ж�������·��
    int k , i , j; 
 
    for( i = 0 ; i <= numnode ; ++i ){
        for( j = 0 ; j <= numnode ; ++j ){
            if(dis[i][j] < INF)
                pathmatirx[i][j] = j;               //��ʼ��·������
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
                    //��������±�Ϊk����·����ԭ�����·������
                    dis[i][j] = dis[i][k] + dis[k][j];   //���µ�ǰ�����Ȩֵ
                    pathmatirx[i][j] = pathmatirx[i][k]; //·������Ϊ�����±�ΪK�Ķ���
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
 
void gen_quicksort(gen *number,int left,int right)//�����������ڽṹ������
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
            //pathmatirx[i][j] = j;               //��ʼ��·������
            if(i == j ) dis[i][j] = 0;
            else dis[i][j] = INF;
        }
    }
    for(i = 0 ; i < MAXN ; ++i ){
        ismustedge[i].flag = false;
    }
}
 
int randsign(float p)//������p����1 
{ 
    if(rand() > (p * 32768)) 
        return 0; 
    else return 1; 
} 
int randbit(int i, int j)//������i��j֮���һ������� 
{ 
    int a , l; 
    l = j - i + 1; 
    a = i + rand() * l / 32768; 
    return a; 
}
