function  Write2TXT(mode, fname, string, value )

fp   =  fopen( fname, mode );


if exist('string')
    
    fprintf( fp, '%s\n\n', string );
    
end

if exist( 'value' )    
    fprintf( fp, '%6.2f       %6.4f\n', value );
end

fclose(fp);