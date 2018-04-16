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
			if(c[i][j] != -1){		//�ߴ��� 
				d[i][j] = c[i][j];
				prev[i][j] = i;
			}	
			else{
				d[i][j] = 9999;	
				prev[i][j] = 0;		//ǰ����㲻���� 
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
		int d=0;	//��¼·������
		
		while(s!=v&&prev[v][u]!=0){				
			way[d++] = prev[v][s];	
			s = prev[v][s];	//ǰ����� 	
		} 		
		for(int i=d-1;i>=0;i--){
			cout<<way[i]-1<<" ->";
			count++;
		}
		return count;				
}

int main()  
{  
	//��Ҫ����Ĳ��������������һ����������������� 
    int n1;								//n1������  
    int k;								//�ع������Ŀ 
	int m; 								//�ع��ߵ���Ŀ  
  	int Danger_num ;						//Σ�ձߵ���Ŀ 
	int Danger_Node[MAX]; 					//Σ�ձ�
  	int c[MAX][MAX];   						//ͼ�����Ӿ���  
  	int mustNode[MAX];						//���бع��� 
  	int most_Node ; 						//���ֻ�ܾ����ĵ� 
   
   	//���� 
    int n2=k+2*m+2;							//���бع������Ŀ���ع���+�ع��ߵ����ˣ�	��������ֹ�㡿 
	int g[MAX][MAX],x[MAX],bestx[MAX];		//��ͼ�����Ӿ��󣬵�ǰ·��������·�� 
	int cl=0,bestl=MAX;						//��ǰ���ѣ����ٻ��� 
    int d[MAX][MAX];						//��������֮������·�� 
    int prev[MAX][MAX]; 					//���ٻ���-ǰ����� 
    int prev2[MAX][MAX];					//���ٽ��-ǰ����� 
	int path[MAX];  							//�������·��	
  	
  	int count1=0;
  	int count2=0;
  	
  	//int d2[MAX][MAX];
   	freopen("����.txt","r",stdin);  
   	//freopen("in160.txt","r",stdin); 
    freopen("���.txt","w",stdout);
    
  	cin>>n1;								//������
	cin>>k;									//�ع����  	
	cin>>m; 								//�ع��� 
	cin>>Danger_num;						//Σ�ձ�   	
	cin>>most_Node;							//���� 	
	
	for(int i=1; i<=n1; i++)  
        for(int j=1; j<=n1; j++)  
            cin>>c[i][j];  	

	mustNode[1]=1;		//����ʱ��Ϊ������㡿 
	for(int i=2;i<=k+1;i++){
	 	cin>>mustNode[i];
		mustNode[i]++;		 				//��Ϊ��ĿΪ0~17���㷨��1~18������ +1 �������ʱ�� -1 
	} 


	n2=k+2*m+2;
	mustNode[n2]=18; 	//���յ㡿 
	for(int i=k+2;i<=n2-1;i+=2){
		cin>>mustNode[i]>>mustNode[i+1];				//1, 2 ~ k+1, k+2 ~ n2-1, n2 
		mustNode[i]++;
		mustNode[i+1]++;			
	}
	

	for(int i=1;i<=2*Danger_num;i+=2) {
	   	cin>>Danger_Node[i]>>Danger_Node[i+1]; 
		c[++Danger_Node[i]][++Danger_Node[i+1]] = NoEdge;	//��ΪNoEdge
			
	}



    for(int i=1; i<=n2; i++){  					//��ʼ��x��bestx 
        x[i]=i;  
        bestx[i]=0;  
    }  

	Floyd(n1,c,prev,d);

	 

	for(int i=1;i<=8;i++){
    	for(int j=1;j<=8;j++){
    		g[i][j]=d[mustNode[i]][mustNode[j]];
		}
	}
	for(int i=k+2;i<=n2-1;i+=2){		//�ع��� ����Ϊ 0 
		g[i][i+1]=g[i+1][i] = 0;
	}

	bestl=TSP(2,n2,cl,bestl,x,bestx,g);
	

	
	for(int i=1;i<=n2;i++){
		path[i] = mustNode[bestx[i]]-1;		
	}			

	int last_path[MAX];
	
	cout<<"���ٻ���·����\n"; 
	for(int i=1;i<n2;i++){
		count1+=out_path(path[i]+1,path[i+1]+1,prev,0);
		
	}
	cout<<path[n2]<<endl;
	
	cout<<"һ��ʹ�õĽ������"<<++count1<<endl;
	
 
    cout<<"��·������Ϊ��";  			//������Ϊ��������ع��ߡ� 
    
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
		
   	for(int i=1; i<=n2; i++){  					//��ʼ��x��bestx 
        x1[i]=i;  
        bestx1[i]=0;  
    }  

	leastNode_Floyd1(n1,node,d1,prev1,c);

	

	for(int i=1;i<=8;i++){
    	for(int j=1;j<=8;j++){
    		g1[i][j]=node[mustNode[i]][mustNode[j]];
		}
	}
	for(int i=k+2;i<=n2-1;i+=2){		//�ع��� ����Ϊ 0 
		g1[i][i+1]=g1[i+1][i] = 0;
	}

	bestl1=TSP(2,n2,cl1,bestl1,x1,bestx1,g1);

	bestl1=0;
	for(int i=1;i<=n2;i++){
		path1[i] = mustNode[bestx1[i]];	
			
	}
	cout<<"\n���ٽ��·����\n";
	for(int i=1;i<n2;i++){
		count2+=out_path(path1[i],path1[i+1],prev1,0);

	}
	cout<<path[n2]<<endl;

	for(int i=1;i<n2;i++){

		bestl1+=d1[path1[i]][path1[i+1]];
	}
	cout<<"һ��ʹ�õĽ������"<<++count2<<endl;
	cout<<"��·������Ϊ��"<<bestl1<<endl;

	cout<<endl; 
	if(count2<=most_Node){
		cout<<"�������Ž�"<<endl; 
		if(count1<=most_Node){
			cout<<"���ٻ���·����Ϊ����·��\n";
		}
		else{
		
			cout<<"����·��Ϊ���ٽ��·��"; 
		}
	}
	else{
		cout<<"�޽�\n";
		if(bestl<bestl1)
		cout<<"����·�������ٻ���·��\n"; 
		else
		cout<<"����·�������ٽ��·��\n";
	} 
	//system("pause"); 
    return 0;  
}  


