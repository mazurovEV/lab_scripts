
# coding: utf-8

# In[150]:


import xml.etree.ElementTree as ET
import random


# In[227]:


class PoAPowersGenerator:
    
    LINE_LENGTH = 60
    
    def __init__(self, xml):
        tree = ET.parse(xml)
        self.root = tree.getroot()
        
        self.current_line = 0
        self.current_end_of_prev_token = 0
    
    def find_all_frames(self, frame=None, offset='', markup=[]):
        #if frame.attrib['name'] != 'Блоки':  # Не ясно почему у фрейма 'Блоки' нет limitations
        #    print(offset + frame.attrib['name'] + 
        #              ", limitations: minOccurs - " + frame.attrib['minOccurs'] + 
        #              ", maxOccurs - " + frame.attrib['maxOccurs'])
        #else:
        #    print(offset + frame.attrib['name'])

        # Тут для каждого фрейма надо понять сколько их будет и делать цикл
        
        if not frame:
            frame = self.root
        
        frame_count = get_tag_count(frame)
        print(str(frame_count) + " блоков " + frame.attrib['name'])
        for i in range(frame_count):
            nodes = []
            frame_markup = {frame.attrib['name']: nodes}
            markup.append(frame_markup)

            for sub_frame in frame.find('children').findall('frame'):
                find_all_frames(sub_frame, offset + "\t", nodes)

            find_quotes(frame, offset, nodes)
            find_entities(frame, offset, nodes)
            find_categories(frame, offset, nodes)
            
        return markup
        
    # TODO: по данным получить распределения кол-ва фреймов и генерить число из них
    # сейчас можно использовать равномерное распределение, но возможно лучше использвать Пуассона

    #В принципе, для наших целей не важно количество, но вообще все это можно было бы использовать для 
    #аугментации прав
    def get_tag_count(self, frame):
        if frame.attrib['name'] == 'Полномочия':
            max_occurs = frame.attrib['maxOccurs']
            return random.randint(int(frame.attrib['minOccurs']), int(max_occurs) if max_occurs != 'unbounded' else 2)
        if frame.attrib['name'] == 'Действие':
            max_occurs = frame.attrib['maxOccurs']
            return random.randint(int(frame.attrib['minOccurs']), int(max_occurs) if max_occurs != 'unbounded' else 2)

        return 1
    
    def find_quotes(self, elem, offset, nodes):
        for quote in elem.find('children').findall('quote'):
            for i in range(get_tag_count(quote)):
                quote_markup = {"value": quote.attrib['name'], 
                                "entity_id": 0, 
                                "text_pos": get_token_coords(quote, current_line, current_end_of_prev_token)}
                nodes.append({quote.attrib['name']: quote_markup})

            #print(offset + "\t" + "quote: " + quote.attrib['name'] + 
            #      ", limitations: minOccurs - " + quote.attrib['minOccurs'] + 
            #      ", maxOccurs - " + quote.attrib['maxOccurs'])

    def get_token_coords(self, tag):
        token_length = len(tag.attrib['name'])

        start_pos = [self.current_line, self.current_end_of_prev_token]

        possible_end_pos = current_end_of_prev_token + token_length
        if possible_end_pos > LINE_LENGTH:
            self.current_line = current_line + 1
            self.current_end_of_prev_token = possible_end_pos - LINE_LENGTH

        end_pos = [self.current_line, self.current_end_of_prev_token]

        return [start_pos, end_pos]
    
    def find_entities(self, elem, offset, nodes):
        for entity in elem.find('children').findall('entity'):
            for i in range(get_tag_count(entity)):
                entity_markup = {"value": "value", "entity_id": 0, "text_pos": [[], []]}
                nodes.append({entity.attrib['name']: entity_markup})

            #print(offset + "\t" + "entity: " + entity.attrib['name'] + 
            #      ", limitations: minOccurs - " + entity.attrib['minOccurs'] + 
            #      ", maxOccurs - " + entity.attrib['maxOccurs'])
            
    def find_categories(self, elem, offset, nodes):
        for category in elem.find('children').findall('category'):
            for i in range(get_tag_count(category)):
                nodes.append({category.attrib['name'] + "__category_list": "category"})

            #print(offset + "\t" + "category: " + category.attrib['name'] + 
            #      ", limitations: minOccurs - " + category.attrib['minOccurs'] + 
            #      ", maxOccurs - " + category.attrib['maxOccurs'])


# In[228]:


generator = PoAPowersGenerator('/home/mazurovev/powers.xml')


# In[229]:


markup = generator.find_all_frames()


# In[196]:


all_markup

