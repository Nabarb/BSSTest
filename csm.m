function c=csm(clean_sig,bad_sig)
for i=1:size(bad_sig,1)
c(1,i)=abs((dot(bad_sig(i,:),clean_sig(i,:))/(norm(clean_sig(i,:))*(norm(bad_sig(i,:))))));
end