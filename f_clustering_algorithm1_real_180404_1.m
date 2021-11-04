function [result1 result2] = f_clustering_algorithm1_real_180404_1(CM,n_m,n_m_g)


    
    mean_CM = mean(mean(CM));

    CM_outlier = sum(CM) - 1;

    mean_CM_out = mean(CM_outlier);

    sd = sqrt(var(CM_outlier));

    CM_out_nd = (CM_outlier-mean_CM_out)/sd;


    [a b] = find(CM_out_nd > 3 | CM_out_nd < -3);

    outlier_i = n_m(b);
    outlier_g = n_m_g(outlier_i);
    % 
    n_m_t = n_m;
    n_m_g_t = n_m_g(n_m_t);
    
    

    n_m_g_t(b) = [];
    n_m_t(b) = [];
    
    
   


result1{1} = n_m_t;
result1{2} = outlier_i;

result2{1} = n_m_g_t;
result2{2} = outlier_g;



end


