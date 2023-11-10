# Testing of cocopaths



## Testing AFPs were resulting domain seq should fold like the AFP



### 1. AFP 

```bash
.
()
.()
(())
.()()*
(()())


Result:  m0*  l0	a m0 b l1	b* m0* a* l2	 c m0 d l3	 d* m0* c* l4	  m0  l5
```



### 2. AFP 
```bash
.
()
().
(())


Result: a* b* m0* c* d* l0	c m0 b l1	 m0*  l2	d c m0 b a l3

```





### 3. AFP 

```bash
.
()
().
()()
().()
()(())



Result:  m0*  l0	 m0  l1	m1*  l2		a m1 b l3	b* m1* a* l4	 m1  l5
```

Like #2 but with first to nodes are added



### 4. AFP

```bash
.
()
().
()()
()().
()()()
()()().

Result:  m0*  l0	m0  l1		m1*  l2		m1  l3		m2*  l4	 m2  l5		m3*  l6
```



Test for seeing if seperated graphs get their own middle domain



### 5. AFP 

```bash
.
()
.()
(())
(()).
(()())


Result: a* b* m0* c* d* l0	e m0 f l1	f* m0* e* l2	c m0 b l3	 m0*  l4	d c m0 b a l5
```


### 6. AFP

AFP where each node gets refolded each step


```bash
.
()
.()
()()
.()()

Result:  m0*  l0	a m0 b l1	c* b* m0* a* d* l2	e d a m0 b c f l3	f* c* b* m0* a* d* e* l4
```


