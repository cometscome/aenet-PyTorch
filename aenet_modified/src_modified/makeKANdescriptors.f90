program main
    use makeKAN,only:initialize,makeKANdescriptor
    implicit none
    character(len=1024) :: infile,outfile
    integer::npoints
    logical::frombinary

    call initialize(infile,outfile,npoints,frombinary) 
    call makeKANdescriptor(infile,outfile,npoints,frombinary)

end program