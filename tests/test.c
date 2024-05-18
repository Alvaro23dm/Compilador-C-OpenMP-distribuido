#pragma omp cluster

#pragma omp target

#pragma omp targgget

#pragma omp target map (tofrom : x) 

#pragma omp cluster broad (width , height)

#pragma omp cluster gather ( rgb[height*width] : chunk(width) )

#pragma omp cluster scatter(matA[nfilas][ncols]:chunk(n))

#pragma omp cluster broad( matB[fB][cB] )

#pragma omp cluster alloc ( width, height )

#pragma omp cluster allgather (matB[nfil][ncols]:chunk(n))

#pragma omp cluster halo (width : chunk(n))

#pragma omp cluster reduction ( + : width )

#pragma omp cluster allreduction ( * : height )

#pragma omp declare cluster

#pragma omp end declare cluster

#pragma omp cluster alloc (width, height, depth)

#pragma omp cluster data broad (width , height)

#pragma omp cluster data gather ( rgb[height*width] : chunk(width) )

#pragma omp cluster data scatter(matA[nfilas][ncols]:chunk(n))

#pragma omp cluster data broad( matB[fB][cB], load )

#pragma omp cluster data allgather (matB[nfil][ncols]:chunk(n))

#pragma omp cluster data halo (width : chunk(n))

#pragma omp cluster data reduction ( || : width )

#pragma omp cluster data allreduction ( * : height )

#pragma omp cluster update broad (height)

#pragma omp cluster update gather ( rgb[height*width] : chunk(width) )

#pragma omp cluster update scatter(matA[nfilas][ncols]:chunk(n))

#pragma omp cluster update broad( matB[fB][cB] )

#pragma omp cluster update allgather (pos[np*nd]:chunk(nd))

#pragma omp cluster update broad (c_max)

#pragma omp cluster update allgather (matB[nfil][ncols]:chunk(n))

#pragma omp cluster update halo (width : chunk(n))

#pragma omp cluster update reduction ( - : width )

#pragma omp cluster update allreduction ( * : height )

#pragma omp cluster teams num_teams (n)  

#pragma omp cluster teams num_teams (n) private (x,y,z) 

#pragma omp cluster teams firstprivate (i,j,k) reduction (-: sum)

#pragma omp cluster teams shared (n,m) thread_limit (15)

#pragma omp cluster distribute collapse (n)

#pragma omp cluster distribute firstprivate(x,y,z) private (a,b,c)

#pragma omp cluster distribute lastprivate (n,m)

#pragma omp cluster teams distribute num_teams (n)  

#pragma omp cluster teams distribute num_teams (n) private (x,y,z) 

#pragma omp cluster teams distribute firstprivate (i,j,k) reduction (-: sum)

#pragma omp cluster teams distribute shared (n,m) thread_limit (15)

#pragma omp cluster teams distribute collapse (n)

#pragma omp cluster teams distribute firstprivate(x,y,z) private (a,b,c)

#pragma omp cluster teams distribute lastprivate (n,m)

#pragma omp cluster teams master

#pragma omp task_async depend (in: x, v[i], w[10:10])

#pragma omp task_async depend (in:c[i:1]) depend(out:d[i:1])

#pragma omp task_async depend (in: x, n[i:3]) depend(out: v[i],w[10:10])


