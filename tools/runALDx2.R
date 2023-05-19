DA_aldx2 = function(phylo, mm_form, effectColNumber = 1, TestColNumber=18, cutoff=0.05, taxlevel="Genus"){
  
  meta = data.frame(sample_data(phylo))
  mm = model.matrix(formula(mm_form), meta)
  
  taxtab = tax_table(phylo)
  rownames(phylo@otu_table@.Data) = taxtab[, taxlevel]
  # input phyloseq, model
  x = aldex.clr(data.frame(otu_table(phylo)), mm, denom = "all")
  glm.text = aldex.glm(x, mm)
  glm.effect = aldex.glm.effect(x)

  print("Test names")
  print(colnames(glm.text))
  print("Effect names")
  print(names(glm.effect))
  
  sig.val = glm.text[, TestColNumber] < cutoff
  sig.val.effect= glm.effect[[effectColNumber]]$effect[sig.val]
  name.effect = names(glm.effect)[effectColNumber]
  
  print(sum(sig.val))
  if(sum(sig.val)==0){
    print("No significant taxa detected")
    return(0)
  } else {
    df = data.frame(effect_size=sig.val.effect, taxa = rownames(glm.effect[[effectColNumber]])[sig.val], var = name.effect)
    
    return(df)
  }
  
  
}