function h_value = Hash(Data, Length, seed)

   Opt = struct('Method', 'MD2', 'Format', 'double');
   Data = [Data , seed];
   temp = DataHash(Data, Opt);
   h_value = mod(sum(temp), Length) + 1; 
end