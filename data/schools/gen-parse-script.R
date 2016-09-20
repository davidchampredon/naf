ncpu <- 4
s.start <- 1
s.end <- 200

batchsize <- floor((s.end-s.start)/ncpu)
a<-vector()
b<-vector()
tmp1 <- s.start
tmp2 <- s.start+batchsize
i <-1
while (tmp2<=s.end) {
	a[i] <- tmp1
	b[i] <- tmp2
	tmp1<-tmp2+1
	tmp2 <- tmp2+batchsize+1
	i<-i+1
}

base.cmd <- 'Rscript hack-website.R'

z <- paste(base.cmd,a,b,'&')
write(z,file='parse_url')
system('chmod +x parse_url')
