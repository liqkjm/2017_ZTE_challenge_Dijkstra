#include<iostream>  
using namespace std;  

#define NoEdge -1 
const int MAX= 200;

void output(int m,int n,int array[MAX]);
 
int TSP(int t,int n,int cl,int bestl,int x[MAX],int bestx[MAX],int g[MAX][MAX] ){  
    int j;  
    if(t>n&&x[n]==n)  
    {  
        //output(1,8,x);
        
        if(cl<=bestl){
            for(j=1; j<=n; j++)  
                bestx[j]=x[j];  
            bestl=cl;       	
		}   
        return bestl;
    }  
    else    
    {  
        for(j=t; j<=n; j++)
        {  
            if(g[x[t-1]][x[j]]!=-1 /*&& (cl+g[x[t-1]][x[j]]<=bestl)*/)
            {  
                swap(x[t],x[j]);       
                cl+=g[x[t-1]][x[t]];   
				bestl = TSP(t+1,n,cl,bestl,x,bestx,g);              
				cl-=g[x[t-1]][x[t]];  
                swap(x[t],x[j]);  
            }  
        }  
    } 
	return bestl;
}  

void Floyd(int n,int c[MAX][MAX],int prev[MAX][MAX],int d[MAX][MAX]){
	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++){
			if(c[i][j] != -1){		//边存在 
				d[i][j] = c[i][j];
				prev[i][j] = i;
			}	
			else{
				d[i][j] = 9999;	
				prev[i][j] = 0;		//前驱结点不存在 
			} 
								
			if(i == j)
				d[i][j] = 0;			
		}
	}	
	for(int k=1;k<=n;k++)				
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++)
				if(d[i][j] > (d[i][k]  + d[k][j])){
					d[i][j] = d[i][k] + d[k][j];
					prev[i][j]= prev[k][j];					
				}	
} 

void leastNode_Floyd1(int  n,int node[MAX][MAX],int d[MAX][MAX],int prev[MAX][MAX], int c[MAX][MAX]){

	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++){
			if(c[i][j] !=NoEdge){
				node[i][j] = 2;
				prev[i][j] = i;
				d[i][j] = c[i][j];
			}
			else{
				node[i][j] = 999;
				d[i][j] = 999;
				prev[i][j] = 0;
			}
			if(i==j){
				node[i][j] = 1;
				d[i][j] = 0;
			}
		}
	}
	for(int k=2;k<n;k++)			
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++){
				if(node[i][j] > (node[i][k]  + node[k][j]-1)){
					
					d[i][j] = d[i][k]+d[k][j];
					node[i][j] = node[i][k] + node[k][j] - 1;
					prev[i][j] = prev[k][j];
				}
				else if(node[i][j]==(node[i][k]+node[k][j]-1)){
					if(d[i][j]>(d[i][k]+d[k][j])){
					d[i][j] = d[i][k]+d[k][j];
					node[i][j] = node[i][k] + node[k][j] - 1;
					prev[i][j] = prev[k][j];	
					}
				}					
			}
}

void output(int m,int n,int Array[MAX][MAX]){
	for(int i=m;i<=n;i++){
		for(int j=m;j<=n;j++){
			cout<<Array[i][j]<<" ";
		}
		cout<<endl;
	}
}

void output(int m,int n,int Array[MAX]){
	for(int i=m;i<=n;i++)
	cout<<Array[i]<<" ";
	cout<<endl;
} 

int out_path(int v,int u,int prev[MAX][MAX],int count){
		// v  1 ~ n 
		//count=0;
		if(prev[v][u]==0)
			cout<<"no path\n";

		int way[18];
		int s=u;	//start
		int d=0;	//记录路径长度
		
		while(s!=v&&prev[v][u]!=0){				
			way[d++] = prev[v][s];	
			s = prev[v][s];	//前驱结点 	
		} 		
		for(int i=d-1;i>=0;i--){
			cout<<way[i]-1<<" ->";
			count++;
		}
		return count;				
}

int main()  
{  
	//需要输入的参数，能清楚描述一个问题的所有条件。 
    int n1;								//n1个顶点  
    int k;								//必过点的数目 
	int m; 								//必过边的数目  
  	int Danger_num ;						//危险边的数目 
	int Danger_Node[MAX]; 					//危险边
  	int c[MAX][MAX];   						//图的连接矩阵  
  	int mustNode[MAX];						//所有必过点 
  	int most_Node ; 						//最多只能经过的点 
   
   	//其他 
    int n2=k+2*m+2;							//所有必过点的数目（必过点+必过边的两端）	【还有起止点】 
	int g[MAX][MAX],x[MAX],bestx[MAX];		//简化图的连接矩阵，当前路径，最优路径 
	int cl=0,bestl=MAX;						//当前花费，最少花费 
    int d[MAX][MAX];						//任意两点之间的最短路径 
    int prev[MAX][MAX]; 					//最少花费-前驱结点 
    int prev2[MAX][MAX];					//最少结点-前驱结点 
	int path[MAX];  							//最终输出路径	
  	
  	int count1=0;
  	int count2=0;
  	
  	//int d2[MAX][MAX];
   	freopen("输入.txt","r",stdin);  
   	//freopen("in160.txt","r",stdin); 
    freopen("输出.txt","w",stdout);
    
  	cin>>n1;								//顶点数
	cin>>k;									//必过点点  	
	cin>>m; 								//必过边 
	cin>>Danger_num;						//危险边   	
	cin>>most_Node;							//最多点 	
	
	for(int i=1; i<=n1; i++)  
        for(int j=1; j<=n1; j++)  
            cin>>c[i][j];  	

	mustNode[1]=1;		//【暂时人为设置起点】 
	for(int i=2;i<=k+1;i++){
	 	cin>>mustNode[i];
		mustNode[i]++;		 				//因为题目为0~17，算法是1~18，所以 +1 ，输出的时候 -1 
	} 


	n2=k+2*m+2;
	mustNode[n2]=18; 	//【终点】 
	for(int i=k+2;i<=n2-1;i+=2){
		cin>>mustNode[i]>>mustNode[i+1];				//1, 2 ~ k+1, k+2 ~ n2-1, n2 
		mustNode[i]++;
		mustNode[i+1]++;			
	}
	

	for(int i=1;i<=2*Danger_num;i+=2) {
	   	cin>>Danger_Node[i]>>Danger_Node[i+1]; 
		c[++Danger_Node[i]][++Danger_Node[i+1]] = NoEdge;	//设为NoEdge
			
	}



    for(int i=1; i<=n2; i++){  					//初始化x，bestx 
        x[i]=i;  
        bestx[i]=0;  
    }  

	Floyd(n1,c,prev,d);

	 

	for(int i=1;i<=8;i++){
    	for(int j=1;j<=8;j++){
    		g[i][j]=d[mustNode[i]][mustNode[j]];
		}
	}
	for(int i=k+2;i<=n2-1;i+=2){		//必过边 设置为 0 
		g[i][i+1]=g[i+1][i] = 0;
	}

	bestl=TSP(2,n2,cl,bestl,x,bestx,g);
	

	
	for(int i=1;i<=n2;i++){
		path[i] = mustNode[bestx[i]]-1;		
	}			

	int last_path[MAX];
	
	cout<<"最少花费路径：\n"; 
	for(int i=1;i<n2;i++){
		count1+=out_path(path[i]+1,path[i+1]+1,prev,0);
		
	}
	cout<<path[n2]<<endl;
	
	cout<<"一共使用的结点数："<<++count1<<endl;
	
 
    cout<<"其路径长度为：";  			//加上置为零的两条必过边。 
    
	int sum=0;
	for(int i=k+2;i<=n2-1;i+=2){
    	sum+=c[mustNode[i]][mustNode[i+1]];
	}
	 
    cout<<bestl+sum<<endl;  

  	int x1[MAX],bestx1[MAX];
  	int d1[MAX][MAX];
  	int g1[MAX][MAX];
  	int prev1[MAX][MAX];
	int cl1=0;
	int bestl1=MAX;
	int node[MAX][MAX];
	int path1[MAX];
		
   	for(int i=1; i<=n2; i++){  					//初始化x，bestx 
        x1[i]=i;  
        bestx1[i]=0;  
    }  

	leastNode_Floyd1(n1,node,d1,prev1,c);

	

	for(int i=1;i<=8;i++){
    	for(int j=1;j<=8;j++){
    		g1[i][j]=node[mustNode[i]][mustNode[j]];
		}
	}
	for(int i=k+2;i<=n2-1;i+=2){		//必过边 设置为 0 
		g1[i][i+1]=g1[i+1][i] = 0;
	}

	bestl1=TSP(2,n2,cl1,bestl1,x1,bestx1,g1);

	bestl1=0;
	for(int i=1;i<=n2;i++){
		path1[i] = mustNode[bestx1[i]];	
			
	}
	cout<<"\n最少结点路径：\n";
	for(int i=1;i<n2;i++){
		count2+=out_path(path1[i],path1[i+1],prev1,0);

	}
	cout<<path[n2]<<endl;

	for(int i=1;i<n2;i++){

		bestl1+=d1[path1[i]][path1[i+1]];
	}
	cout<<"一共使用的结点数："<<++count2<<endl;
	cout<<"其路径长度为："<<bestl1<<endl;

	cout<<endl; 
	if(count2<=most_Node){
		cout<<"存在最优解"<<endl; 
		if(count1<=most_Node){
			cout<<"最少花费路径即为最优路径\n";
		}
		else{
		
			cout<<"最优路径为最少结点路径"; 
		}
	}
	else{
		cout<<"无解\n";
		if(bestl<bestl1)
		cout<<"次优路径：最少花费路径\n"; 
		else
		cout<<"次有路径：最少结点路径\n";
	} 
	//system("pause"); 
    return 0;  
}  


