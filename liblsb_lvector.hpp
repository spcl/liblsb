/*
* Copyright (c) 2015 ETH-Zurich. All rights reserved.
* Use of this source code is governed by a BSD-style license that can be
* found in the LICENSE file.
*/

#ifndef __LVECTOR_HPP__
#define __LVECTOR_HPP__

//#include <stdexcept>
#include <vector>
#include <stdio.h>

template <typename T>
struct lv_entry{
    std::vector<T> vect;
    lv_entry * next;
    //unsigned int id;
    unsigned int l,r;
    lv_entry(unsigned int l, size_t vsize){
        this->next = NULL;
        this->l = l;
        this->r = l + vsize - 1;
        vect.reserve(vsize);
    }

    int isAtCapacity(){
        return vect.capacity()<=vect.size();
    }
};


template <typename T>
class lvector{

private:
    lv_entry<T> * list; // list head pointer
    lv_entry<T> * list_tail; // list tail pointer
    lv_entry<T> * last_pointer; // caches the last pointed entry for reading
    lv_entry<T> * list_insert; //the block where insert new data
    size_t _size; // actual size of the vector
    size_t _capacity; // block id counter
public:

    lvector(){
        this->_size=0;
        this->_capacity = 1;
        list = new lv_entry<T>(0, _capacity);
        list_tail=list;
        last_pointer=list;
        list_insert = list;
        list->next=NULL;
    }

    ~lvector(){
        while (list!=NULL){
            lv_entry<T> * tmp = list->next;
            delete list;
            list = tmp;
        }
    }
    

    inline void extend(size_t toadd){
        lv_entry<T> * newentry = new lv_entry<T>(_capacity, toadd);
        _capacity += toadd;
        list_tail->next = newentry;
        list_tail = newentry;
    }

    void clear(){
        lv_entry<T> *ptr, *tmp;
        ptr = list->next;
        while (ptr!=NULL){
            tmp = ptr->next;
            delete ptr;
            ptr = tmp;
        }

        int oldcap = list->vect.capacity();
        list->vect.clear();
        list->next = NULL;
        _capacity=1;
        _size=0;
        list_insert=list;
        list_tail=list;
        last_pointer=list;
    }

    T& back(){
        return list_insert->vect.back();
    }

    void push_back(T el){

        //this is ugly but I need it since someone could have extended the list with reserve()
        while (list_insert!=NULL && list_insert->isAtCapacity()) {
            //printf("lsb_lvector: advancing in push_back\n");
            list_insert = list_insert->next;
        }

        if (list_insert==NULL){
            extend(_capacity);
            list_insert = list_tail;
        }

        list_insert->vect.push_back(el);
        this->_size++;
    }

    void reserve(size_t s){
        if (s>_capacity) extend(s - _capacity);     
    }

    //int i=0;
    T& operator[](unsigned int idx){
        //unsigned int rid = (unsigned int) (idx / (float) vsize);
        //printf("lib: idx: %i; rid: %i;\n", idx, rid);
        //if (last_pointer->id == rid) printf("reusing last pointer for %i\n", idx);
        //int c = 0;
        if (last_pointer->l > idx) last_pointer = list;
        while (last_pointer!=NULL && (last_pointer->r < idx)) {
            last_pointer = last_pointer->next;
            //c++;
        }

        /*if (c>0){
            printf("Walked %i steps after %i\n", c, i);
            i=0;
        }else i++;*/

        //if (last_pointer == NULL) throw std::out_of_range("out_of_range");
        
        //printf("accessing block %i\n----------------\n", last_pointer->id);
        return last_pointer->vect[idx - last_pointer->l];
    }

    size_t size(){ return _size; }

    size_t capacity() { return _capacity; }    
/*
    int check(){
        printf("total capacity: %i\n", (int) this->capacity());
        lv_entry<T> * tmp = list;
        while (tmp!=NULL){
            printf("vect.capacity(): %i;\n", (int) tmp->vect.capacity());
            tmp = tmp->next;
        }
    }
*/
};
 
#endif /* __LVECTOR_HPP__ */
