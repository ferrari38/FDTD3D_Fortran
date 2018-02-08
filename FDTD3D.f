!-----------------------------------------------------------------------
!        FDTDの共通変数
!-----------------------------------------------------------------------
      module fdtd_variable
!解析領域                                
      integer,parameter::nxx=200,nyy=200,nzz=200   !解析領域分割数
      integer,parameter::nstep=300                !計算ステップ数
      real,parameter::dx=0.005,dy=0.005,dz=0.005   !セルサイズ
      real::dt,t                                   !時間ステップ,時間
!全計算領域
      integer,parameter::nx=nxx,ny=nyy,nz=nzz
!電界，磁界の配列，係数の配列
      real::ex(0:nx,0:ny,0:nz),ey(0:nx,0:ny,0:nz),ez(0:nx,0:ny,0:nz)
      real::hx(0:nx,0:ny,0:nz),hy(0:nx,0:ny,0:nz),hz(0:nx,0:ny,0:nz)
      real::aex(0:nx,0:ny,0:nz),aey(0:nx,0:ny,0:nz),aez(0:nx,0:ny,0:nz)
      real::bexy(0:nx,0:ny,0:nz),bexz(0:nx,0:ny,0:nz)
      real::beyx(0:nx,0:ny,0:nz),beyz(0:nx,0:ny,0:nz)
      real::bezx(0:nx,0:ny,0:nz),bezy(0:nx,0:ny,0:nz)
      real::amx(0:nx,0:ny,0:nz),amy(0:nx,0:ny,0:nz),amz(0:nx,0:ny,0:nz)
      real::bmxy(0:nx,0:ny,0:nz),bmxz(0:nx,0:ny,0:nz)
      real::bmyx(0:nx,0:ny,0:nz),bmyz(0:nx,0:ny,0:nz)
      real::bmzx(0:nx,0:ny,0:nz),bmzy(0:nx,0:ny,0:nz)
!媒質定数の配列と背景媒質定数
      real::epsd(-1:nx,-1:ny,-1:nz),sgmed(-1:nx,-1:ny,-1:nz)
      real::mud(-1:nx,-1:ny,-1:nz), sgmmd(-1:nx,-1:ny,-1:nz)
      real,parameter::epsbk=1.0,mubk=1.0,sigebk=0.0,sigmbk=0.0 !背景媒質
!散乱体
      integer,parameter::ic=nx/2, jc=ny/2, kc=nz/2    !散乱直方体の中心
      integer,parameter::lx2=20, ly2=20, lz2=20       !(直方体の寸法）/2
      real,parameter::epsr=5.0                        !直方体の比誘電率
!励振パルス
      real,parameter::duration=0.1e-9,t0=4.0*duration   !パルス幅，ピーク時刻 
      integer,parameter::ifed=ic-lx2-20, jfed=jc, kfed=kc   !給電位置
!励振の種類
      integer::lfeed,lsoft,lhard
      parameter(lsoft=1,lhard=2)
      parameter(lfeed=lhard)
      real::befed                                       !ソフト給電の係数
      real,parameter::dl=0.001                          !ハード給電の間隙
!定数
      real,parameter::eps0=8.854188e-12,mu0=1.256637e-6 !真空の誘電率，透磁率
      real,parameter::c=2.9979246e8                     !光速
      end module fdtd_variable
!-----------------------------------------------------------------------
!       メインプログラム
!-----------------------------------------------------------------------      
      program fdtd_3d
      use fdtd_variable

      call setup()                    !FDTDの初期設定

      t=dt
      do n=1,nstep                    !繰り返し計算
         write(*,*)'Time step:',n
         call e_cal                   !電界の計算
         call feed()                  !電流源の励振
         t=t+0.5*dt                   !時間の更新
         call h_cal                   !磁界の計算
         t=t+0.5*dt                   !時間の更新
         call out_emf(n)              !計算結果の出力
      end do
      end program fdtd_3d
!-----------------------------------------------------------------------
!       計算結果の出力
!-----------------------------------------------------------------------
      subroutine out_emf(n)
      use fdtd_variable
      integer,parameter::io=ic-lx2-40, jo=ny/2, ko=nz/2    !観測点

!出力空間の設定
      is=lpml
      ie=nx-lpml
      js=lpml
      je=ny-lpml
      k=nz/2

      if(n == 1) open(02,file='eztm.txt')
      write(02,111)t,ez(io,jo,ko)                   !観測点の過渡電界

      if(n==160)then
         open(03,file="ez160h.txt")
         do j=js,je
            do i=is,ie
               write(03,222)i-lpml,j-lpml,ez(i,j,k) !電界の空間分布
            end do
         end do
         close(03)
      end if

      if(n == nstep) close(02)
  111 format(2e18.9)
  222 format(2i5,e15.6)

      end subroutine out_emf 

!-----------------------------------------------------------------------
!       励振波源
!-----------------------------------------------------------------------      
      subroutine feed()
      use fdtd_variable
      real::iz


      ez(ifed,jfed,kfed)=exp(-((t-t0)/duration)**2)/dl   !ハード給電

      end subroutine feed
!-----------------------------------------------------------------------
!     FDTDの初期設定：
!-----------------------------------------------------------------------       
      subroutine setup()
      use fdtd_variable
      real::mu
      integer::i,j,k

!時間ステップ
      dt=0.99999/(c*sqrt(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz))) 
!背景媒質
      do k=-1,nz
         do j=-1,ny
            do i=-1,nx
               epsd(i,j,k)=epsbk
               mud(i,j,k)=mubk
               sgmed(i,j,k)=sigebk
               sgmmd(i,j,k)=sigmbk
            end do
         end do
      end do

!散乱直方体の媒質定数（背景媒質に上書き）
      call epsmu()

!波源の位置における係数
      eps=0.25*(epsd(ifed,jfed,kfed)+epsd(ifed-1,jfed,kfed)         
     &          +epsd(ifed,jfed-1,kfed)+epsd(ifed-1,jfed-1,kfed))*eps0
      befed=dt/eps

!係数の計算
      do k=0,nz
        do j=0,ny
          do i=0,nx
            eps=0.25*(epsd(i,j,k)+epsd(i,j-1,k)
     &            +epsd(i,j,k-1)+epsd(i,j-1,k-1))*eps0 
            sgm=0.25*(sgmed(i,j,k)+sgmed(i,j-1,k)
     &            +sgmed(i,j,k-1)+sgmed(i,j-1,k-1))
            a=0.5*sgm*dt/eps 
            aex(i,j,k)=(1.0-a)/(1.0+a)
            bexy(i,j,k)=dt/eps/(1.0+a)/dy
            bexz(i,j,k)=dt/eps/(1.0+a)/dz

            eps=0.25*(epsd(i,j,k)+epsd(i-1,j,k)
     &            +epsd(i,j,k-1)+epsd(i-1,j,k-1))*eps0
            sgm=0.25*(sgmed(i,j,k)+sgmed(i-1,j,k)
     &            +sgmed(i,j,k-1)+sgmed(i-1,j,k-1))
            a=0.5*sgm*dt/eps 
            aey(i,j,k)=(1.0-a)/(1.0+a)
            beyx(i,j,k)=dt/eps/(1.0+a)/dx
            beyz(i,j,k)=dt/eps/(1.0+a)/dz

            eps=0.25*(epsd(i,j,k)+epsd(i-1,j,k)
     &            +epsd(i,j-1,k)+epsd(i-1,j-1,k))*eps0
            sgm=0.25*(sgmed(i,j,k)+sgmed(i-1,j,k)
     &            +sgmed(i,j-1,k)+sgmed(i-1,j-1,k))
            a=0.5*sgm*dt/eps 
            aez(i,j,k)=(1.0-a)/(1.0+a)
            bezx(i,j,k)=dt/eps/(1.0+a)/dx
            bezy(i,j,k)=dt/eps/(1.0+a)/dy

            mu=0.5*(mud(i,j,k)+mud(i-1,j,k))*mu0
            sgmm=0.5*(sgmmd(i,j,k)+sgmmd(i-1,j,k))
            a=0.5*sgmm*dt/mu 
            amx(i,j,k)=(1.0-a)/(1.0+a)
            bmxy(i,j,k)=dt/mu/(1.0+a)/dy
            bmxz(i,j,k)=dt/mu/(1.0+a)/dz

            mu=0.5*(mud(i,j,k)+mud(i,j-1,k))*mu0
            sgmm=0.5*(sgmmd(i,j,k)+sgmmd(i,j-1,k))
            a=0.5*sgmm*dt/mu 
            amy(i,j,k)=(1.0-a)/(1.0+a)
            bmyz(i,j,k)=dt/mu/(1.0+a)/dz
            bmyx(i,j,k)=dt/mu/(1.0+a)/dx

            mu=0.5*(mud(i,j,k)+mud(i,j,k-1))*mu0
            sgmm=0.5*(sgmmd(i,j,k)+sgmmd(i,j,k-1))
            a=0.5*sgmm*dt/mu 
            amz(i,j,k)=(1.0-a)/(1.0+a)
            bmzx(i,j,k)=dt/mu/(1.0+a)/dx
            bmzy(i,j,k)=dt/mu/(1.0+a)/dy
          end do
        end do
      end do
      end subroutine setup
      
!-----------------------------------------------------------------------
!　　直方体の媒質定数
!-----------------------------------------------------------------------
      subroutine epsmu()   
      use fdtd_variable 

      do k=kc-lz2, kc+lz2-1
         do j=jc-ly2, jc+ly2-1
            do i=ic-lx2, ic+lx2-1
               epsd(i,j,k)=epsr
               mud(i,j,k)=1.0
               sgmed(i,j,k)=0.0
               sgmmd(i,j,k)=0.0
            end do
         end do
      end do
      end subroutine epsmu

!-----------------------------------------------------------------------
!     電界の計算
!-----------------------------------------------------------------------
      subroutine e_cal()
      use fdtd_variable
      integer::i,j,k
!Ex
      do k=1,nz-1 
         do j=1,ny-1    
            do i=0,nx-1       
               ex(i,j,k)=aex(i,j,k)*ex(i,j,k)+bexy(i,j,k)*(hz(i,j,k)
     &             -hz(i,j-1,k))-bexz(i,j,k)*(hy(i,j,k)-hy(i,j,k-1))
            end do
         end do
      end do
!Ey
      do k=1,nz-1
         do j=0,ny-1
            do i=1,nx-1
               ey(i,j,k)=aey(i,j,k)*ey(i,j,k)+beyz(i,j,k)*(hx(i,j,k)
     &             -hx(i,j,k-1))-beyx(i,j,k)*(hz(i,j,k)-hz(i-1,j,k))
            end do
         end do
      end do
!Ez
      do k=0,nz-1
         do j=1,ny-1
            do i=1,nx-1
               ez(i,j,k)=aez(i,j,k)*ez(i,j,k)+bezx(i,j,k)*(hy(i,j,k)
     &             -hy(i-1,j,k))-bezy(i,j,k)*(hx(i,j,k)-hx(i,j-1,k))
            end do
         end do
      end do
      end subroutine e_cal
!-----------------------------------------------------------------------
!     磁界の計算
!-----------------------------------------------------------------------
      subroutine h_cal()
      use fdtd_variable      
      integer::i,j,k
!Hx
      do k=0,nz-1
         do j=0,ny-1
            do i=1,nx-1
               hx(i,j,k)=amx(i,j,k)*hx(i,j,k)-bmxy(i,j,k)*(ez(i,j+1,k)
     &                   -ez(i,j,k))+bmxz(i,j,k)*(ey(i,j,k+1)-ey(i,j,k))
            end do
         end do
      end do
!Hy
      do k=0,nz-1
         do j=1,ny-1
            do i=0,nx-1
               hy(i,j,k)=amy(i,j,k)*hy(i,j,k)-bmyz(i,j,k)*(ex(i,j,k+1)
     &                   -ex(i,j,k))+bmyx(i,j,k)*(ez(i+1,j,k)-ez(i,j,k))
            end do
         end do
      end do
!Hz
      do k=1,nz-1
         do j=0,ny-1
            do i=0,nx-1
               hz(i,j,k)=amz(i,j,k)*hz(i,j,k)-bmzx(i,j,k)*(ey(i+1,j,k)
     &                   -ey(i,j,k))+bmzy(i,j,k)*(ex(i,j+1,k)-ex(i,j,k))
            end do
         end do
      end do
      end subroutine h_cal
