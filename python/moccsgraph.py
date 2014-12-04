#coding:utf-8

import os, sys, re

import tkgraph
graph = tkgraph.Graph('graph.pdf', tkgraph.Axis(0, 1, 50, 175),
                      tkgraph.Axis(0, 1, 600, 725))
graph.setFont('Courier', 6)

with open(sys.argv[1]) as fi:
    for line in fi:
        items = line.strip().split('\t')
        if len(items) < 4: continue
        graph.drawString(graph.left, graph.top + 10, items[0] + ' ' + items[1])
        path = graph.beginPath()
        count = items[1].split('/')
        total = int(count[0])
        expected = float(count[1])
        values = [float(x_) / total for x_ in items[-1].split(',')]
        path.moveTo(graph.left, graph.bottom)
        acc = 0.0
        for i, value in enumerate(values):
            acc += value
            x = graph.x(float(i + 1) / len(values))
            y = graph.y(acc)
            path.lineTo(x, y)
        #close path
        path.lineTo(graph.right, graph.top)
        path.lineTo(graph.left, graph.bottom)
        path.close()
        graph.setStrokeColor('crimson')
        graph.drawPath(path)
        graph.setStrokeColor('black')
        graph.line(graph.left, graph.bottom, graph.right, graph.top)
        graph.shift()
graph.save()
