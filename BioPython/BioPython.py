#!/usr/bin/env python
# coding: utf-8

# In[85]:


import Bio.Seq
from Bio.Seq import reverse_complement


# In[86]:


dnaseq = "atcggtact"


# In[87]:


seqobj = Bio.Seq.Seq(dnaseq)


# In[88]:


seqobj.translate()


# In[89]:


myseq = "atgcgtagtc"


# In[ ]:





# In[ ]:





# In[146]:


input = 'atg---atatcccgtatatcccgcactgttttacgatcccggatga'


# In[147]:


dna = Bio.Seq.Seq(input)


# In[148]:


output = dna.transcribe() 


# In[149]:


print(output)


# In[150]:


input_rna = 'aug---auaucccguauaucccgcacuguuuuacgaucccggauga'


# In[151]:


rna = Bio.Seq.Seq(input_rna)


# In[152]:


output = rna.back_transcribe()


# In[153]:


print(output)


# In[157]:


input = "atg---atatcccgtatatcccgcactgttttacgatcccggatg"


# In[158]:


dna = Bio.Seq.Seq(input)


# In[159]:


output = dna.translate()


# In[160]:


print(output)


# In[161]:


input = 'atgatatcccgtatat'


# In[162]:


dna = Bio.Seq.Seq(input)


# In[163]:


output = dna.strip("a")


# In[164]:


print(output)


# In[165]:


input = 'aagatatcccgtataaa'


# In[166]:


dna = Bio.Seq.Seq(input)


# In[167]:


output = dna.strip("ag")


# In[168]:


print(output)


# In[169]:


Bio.Seq.Seq('aagatatcccgtataaa').strip("agt")


# In[ ]:




